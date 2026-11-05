function [frames, syncMetric] = decode_audio_stream(rxAudio, cfg, frameKind)

    rxBaseband = audio_to_baseband(rxAudio, cfg);
    [frames, syncMetric] = decode_baseband_stream(rxBaseband, cfg, frameKind);
end

function baseband = audio_to_baseband(audio, cfg)
    audio = double(audio(:));
    t = (0:numel(audio)-1).' / cfg.audioFs;
    mixed = 2 * audio .* exp(-1j * 2 * pi * cfg.audioCarrierHz * t);
    baseband = resample(mixed, cfg.basebandFs, cfg.audioFs);
    baseband = baseband(:);
end

function [frames, syncMetric] = decode_baseband_stream(rxBaseband, cfg, frameKind)
    [frameStartIndices, syncMetric] = detect_preambles(rxBaseband, cfg, frameKind);
    frames = repmat(struct( ...
        'valid', false, ...
        'infoBits', [], ...
        'analysis', []), numel(frameStartIndices), 1);

    for k = 1:numel(frameStartIndices)
        [infoBits, crcError, frameAnalysis] = ...
            decode_phy_frame(rxBaseband, frameStartIndices(k), cfg, frameKind);
        frames(k).valid = crcError == 0;
        frames(k).infoBits = infoBits;
        frames(k).analysis = frameAnalysis;
    end
end

function [frameStartIndices, syncMetric] = detect_preambles(rxBaseband, cfg, frameKind)
    rxBaseband = rxBaseband(:);
    frameParams = frame_params(cfg, frameKind);

    maxStartIndex = numel(rxBaseband) - frameParams.frameBasebandSamples + 1;
    if maxStartIndex < 1
        frameStartIndices = [];
        syncMetric = [];
        return;
    end

    syncMetric = zeros(maxStartIndex, 1);

    for k = 1:maxStartIndex
        firstRepeatedSegment = rxBaseband(k : k + cfg.nFFT - 1);
        secondRepeatedSegment = rxBaseband(k + cfg.nFFT : k + 2 * cfg.nFFT - 1);
        firstEnergy = sum(abs(firstRepeatedSegment).^2);
        secondEnergy = sum(abs(secondRepeatedSegment).^2);
        correlation = sum(conj(firstRepeatedSegment) .* secondRepeatedSegment);
        syncMetric(k) = abs(correlation)^2 / (firstEnergy * secondEnergy + eps);
    end

    frameStartIndices = [];
    remainingMetric = syncMetric;
    peakValue = max([0; syncMetric]);
    candidateThreshold = 0.8 * peakValue;
    minimumFrameSpacing = max(1, floor(0.80 * frameParams.frameBasebandSamples));
    while true
        [localPeak, peakIndex] = max(remainingMetric);
        if localPeak <= 0 || localPeak < candidateThreshold
            break;
        end

        frameStartIndex = peakIndex;
        plateauThreshold = 0.8 * localPeak;
        while frameStartIndex > 1 && syncMetric(frameStartIndex - 1) >= plateauThreshold
            frameStartIndex = frameStartIndex - 1;
        end
        frameStartIndices = [frameStartIndices; frameStartIndex];

        halfSpacing = floor(minimumFrameSpacing / 2);
        suppressedIndices = max(1, peakIndex - halfSpacing) : min(maxStartIndex, peakIndex + halfSpacing);
        remainingMetric(suppressedIndices) = 0;
    end
    frameStartIndices = sort(frameStartIndices);
end

function [infoBits, crcError, analysis] = decode_phy_frame(rxBaseband, frameStartIndex, cfg, frameKind)
    frameParams = frame_params(cfg, frameKind);
    rxBaseband = rxBaseband(:);
    frameEndIndex = frameStartIndex + frameParams.frameBasebandSamples - 1;

    receivedFrame = rxBaseband(frameStartIndex:frameEndIndex);
    knownPreambleSubcarriers = generate_preamble_subcarriers(cfg);
    payloadStartIndex = frameParams.preambleSamples + 1;
    activeGridBeforeCfo = waveform_to_active_grid(receivedFrame, payloadStartIndex, frameParams.nDataSymbols, cfg);

    firstPreambleIndices = cfg.cpLen + (1:cfg.nFFT);
    secondPreambleIndices = cfg.cpLen + cfg.nFFT + (1:cfg.nFFT);
    firstPreambleSymbol = receivedFrame(firstPreambleIndices);
    secondPreambleSymbol = receivedFrame(secondPreambleIndices);
    preambleCorrelation = sum(conj(firstPreambleSymbol) .* secondPreambleSymbol);
    cfoNormalized = angle(preambleCorrelation) / (2 * pi * cfg.nFFT);

    sampleIndices = (0:numel(receivedFrame)-1).';
    cfoCorrectedFrame = receivedFrame .* exp(-1j * 2 * pi * cfoNormalized * sampleIndices);

    firstPreambleSymbol = cfoCorrectedFrame(firstPreambleIndices);
    secondPreambleSymbol = cfoCorrectedFrame(secondPreambleIndices);
    averagedPreambleSymbol = 0.5 * (firstPreambleSymbol + secondPreambleSymbol);
    preambleSpectrum = fftshift(fft(averagedPreambleSymbol, cfg.nFFT) / sqrt(cfg.nFFT));
    receivedPreambleSubcarriers = preambleSpectrum(cfg.activeFftIndices);
    channelEstimate = receivedPreambleSubcarriers ./ knownPreambleSubcarriers;

    activeGridAfterCfo = waveform_to_active_grid(cfoCorrectedFrame, payloadStartIndex, frameParams.nDataSymbols, cfg);
    [compensatedActiveGrid, noiseVariance] = equalize_ofdm_symbols(activeGridAfterCfo, channelEstimate, cfg);

    pilotRatiosBeforeCfo = activeGridBeforeCfo(cfg.continuousPilotIndices, :) ./ cfg.continuousPilotSymbols;
    pilotRatiosAfterCompensation = compensatedActiveGrid(cfg.continuousPilotIndices, :) ./ cfg.continuousPilotSymbols;

    equalizedDataSymbols = reshape(compensatedActiveGrid(cfg.dataActiveIndices, :), [], 1);
    equalizedDataSymbols = equalizedDataSymbols(1:frameParams.encodedBitCount / cfg.bitsPerSymbol);
    llrs = nrSymbolDemodulate(equalizedDataSymbols, cfg.modulation, noiseVariance, 'DecisionType', 'soft');

    [decodedScrambled, crcError] = ldpc_decode_bits(llrs, frameParams.infoBits, frameParams.codeRate, cfg.modulation, frameParams.maxIterations);
    infoBits = scramble_bits(decodedScrambled, frameParams.scramblerSeed);

    analysis.equalizedDataSymbols = equalizedDataSymbols;
    analysis.pilotPhaseBeforeCfo = unwrap(angle(pilotRatiosBeforeCfo), [], 2);
    analysis.pilotPhaseAfterCompensation = unwrap(angle(pilotRatiosAfterCompensation), [], 2);
end

function activeGrid = waveform_to_active_grid(rxSignal, payloadStartIndex, numberOfSymbols, cfg)
    rxSignal = rxSignal(:);
    symbolLength = cfg.nFFT + cfg.cpLen;
    activeGrid = zeros(cfg.nActive, numberOfSymbols);

    for symbolIndex = 1:numberOfSymbols
        cpStartIndex = payloadStartIndex + (symbolIndex-1) * symbolLength;
        fftStartIndex = cpStartIndex + cfg.cpLen;
        fftIndices = fftStartIndex : fftStartIndex + cfg.nFFT - 1;
        usefulSymbol = rxSignal(fftIndices);
        spectrum = fftshift(fft(usefulSymbol, cfg.nFFT) / sqrt(cfg.nFFT));
        activeGrid(:, symbolIndex) = spectrum(cfg.activeFftIndices);
    end
end

function [compensatedActiveGrid, noiseVariance] = equalize_ofdm_symbols(activeGridAfterCfo, channelEstimate, cfg)
    safeChannelEstimate = channelEstimate(:);
    safeChannelEstimate(abs(safeChannelEstimate) < cfg.minimumChannelMagnitude) = 1;
    channelEqualizedGrid = activeGridAfterCfo ./ safeChannelEstimate;

    compensatedActiveGrid = track_continuous_pilots(channelEqualizedGrid, cfg);

    compensatedPilots = compensatedActiveGrid(cfg.continuousPilotIndices, :);
    pilotError = compensatedPilots - cfg.continuousPilotSymbols;
    noiseVariance = max(mean(abs(pilotError(:)).^2), cfg.minNoiseVariance);
end

function compensatedActiveGrid = track_continuous_pilots(channelEqualizedGrid, cfg)
    receivedPilots = channelEqualizedGrid(cfg.continuousPilotIndices, :);
    pilotRatio = receivedPilots ./ cfg.continuousPilotSymbols;
    commonPilotGain = mean(pilotRatio, 1);

    weakPilotSymbols = abs(commonPilotGain) < cfg.minimumChannelMagnitude;
    commonPilotGain(weakPilotSymbols) = 1;

    compensatedActiveGrid = channelEqualizedGrid ./ commonPilotGain;
end

function preambleSubcarriers = generate_preamble_subcarriers(cfg)
    preambleBits = prbs_bits(cfg.nActive, cfg.preambleSeed);
    preambleSubcarriers = cfg.preambleAmplitude * (2 * preambleBits - 1);
end

function [bits, crcError] = ldpc_decode_bits(llrs, transportBlockLength, codeRate, modulation, maxIterations)
    llrs = double(llrs(:));
    cbsInfo = nrDLSCHInfo(transportBlockLength, codeRate);
    rateRecovered = nrRateRecoverLDPC(llrs, transportBlockLength, codeRate, 0, modulation, 1, cbsInfo.C);
    decodedBlocks = nrLDPCDecode(rateRecovered, cbsInfo.BGN, maxIterations);
    [transportBlockWithCrc, ~] = nrCodeBlockDesegmentLDPC(decodedBlocks, cbsInfo.BGN, transportBlockLength + cbsInfo.L);
    [bits, crcError] = nrCRCDecode(transportBlockWithCrc, cbsInfo.CRC);

    bits = double(bits(1:transportBlockLength));
end

function scrambledBits = scramble_bits(bits, seed)
    bits = bits(:) ~= 0;
    seq = prbs_bits(numel(bits), seed) ~= 0;
    scrambledBits = double(xor(bits, seq));
end

function seq = prbs_bits(nBits, seed)
    state = false(15, 1);
    seed = uint32(seed);
    for k = 1:15
        state(k) = bitget(seed, 16 - k) ~= 0;
    end

    seq = false(nBits, 1);
    for k = 1:nBits
        seq(k) = state(end);
        feedback = xor(state(15), state(14));
        state = [feedback; state(1:end-1)];
    end
    seq = double(seq);
end

function frameParams = frame_params(cfg, frameKind)
    switch lower(frameKind)
        case 'data'
            frameParams.infoBits = cfg.dataInfoBits;
            frameParams.codeRate = cfg.dataCodeRate;
            frameParams.scramblerSeed = cfg.scramblerSeedData;
            frameParams.maxIterations = cfg.ldpcMaxIterationsData;
        case 'feedback'
            frameParams.infoBits = cfg.feedbackInfoBits;
            frameParams.codeRate = cfg.feedbackCodeRate;
            frameParams.scramblerSeed = cfg.scramblerSeedFeedback;
            frameParams.maxIterations = cfg.ldpcMaxIterationsFeedback;
    end

    frameParams.codedBitsPerOfdmSymbol = cfg.nDataSubcarriers * cfg.bitsPerSymbol;
    minCodedBits = ceil(frameParams.infoBits / frameParams.codeRate);
    frameParams.nDataSymbols = ceil(minCodedBits / frameParams.codedBitsPerOfdmSymbol);
    frameParams.encodedBitCount = frameParams.nDataSymbols * frameParams.codedBitsPerOfdmSymbol;
    frameParams.preambleSamples = cfg.cpLen + 2 * cfg.nFFT;
    frameParams.symbolSamples = cfg.nFFT + cfg.cpLen;
    frameParams.frameBasebandSamples = frameParams.preambleSamples + frameParams.nDataSymbols * frameParams.symbolSamples;
end
