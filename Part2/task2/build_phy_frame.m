function [frameBaseband, codedBits] = build_phy_frame(infoBits, cfg, frameKind)

    frameParams = frame_params(cfg, frameKind);
    infoBits = double(infoBits(:) ~= 0);

    scrambledBits = scramble_bits(infoBits, frameParams.scramblerSeed);
    codedBits = ldpc_encode_bits(scrambledBits, frameParams.codeRate, cfg.modulation, frameParams.encodedBitCount);
    modulatedSymbols = nrSymbolModulate(codedBits, cfg.modulation);

    activeGrid = zeros(cfg.nActive, frameParams.nDataSymbols);
    for symbolIndex = 1:frameParams.nDataSymbols
        activeSubcarriers = zeros(cfg.nActive, 1);
        activeSubcarriers(cfg.continuousPilotIndices) = cfg.continuousPilotSymbols;
        dataIndices = (symbolIndex - 1) * cfg.nDataSubcarriers + (1:cfg.nDataSubcarriers);
        activeSubcarriers(cfg.dataActiveIndices) = modulatedSymbols(dataIndices);
        activeGrid(:, symbolIndex) = activeSubcarriers;
    end

    preamble = generate_preamble(cfg);
    dataWaveform = active_grid_to_waveform(activeGrid, cfg);
    frameBaseband = [preamble; dataWaveform];
end

function codedBits = ldpc_encode_bits(bits, codeRate, modulation, encodedBitCount)
    bits = double(bits(:) ~= 0);
    cbsInfo = nrDLSCHInfo(numel(bits), codeRate);
    bitsWithCrc = nrCRCEncode(bits, cbsInfo.CRC);
    codeBlocks = nrCodeBlockSegmentLDPC(bitsWithCrc, cbsInfo.BGN);
    encodedBlocks = nrLDPCEncode(codeBlocks, cbsInfo.BGN);
    codedBits = nrRateMatchLDPC(encodedBlocks, double(encodedBitCount), 0, modulation, 1);
    codedBits = double(codedBits(:));
end

function preamble = generate_preamble(cfg)
    preambleBits = prbs_bits(cfg.nActive, cfg.preambleSeed);
    preambleSubcarriers = cfg.preambleAmplitude * (2 * preambleBits - 1);

    shiftedSpectrum = zeros(cfg.nFFT, 1);
    shiftedSpectrum(cfg.activeFftIndices) = preambleSubcarriers;
    preambleSymbol = ifft(ifftshift(shiftedSpectrum), cfg.nFFT) * sqrt(cfg.nFFT);

    preamble = [preambleSymbol(end-cfg.cpLen+1:end); preambleSymbol; preambleSymbol];
end

function waveform = active_grid_to_waveform(activeGrid, cfg)
    nSymbols = size(activeGrid, 2);
    symbolLength = cfg.nFFT + cfg.cpLen;
    waveform = zeros(nSymbols * symbolLength, 1);

    for symbolIndex = 1:nSymbols
        shiftedSpectrum = zeros(cfg.nFFT, 1);
        shiftedSpectrum(cfg.activeFftIndices) = activeGrid(:, symbolIndex);
        usefulSymbol = ifft(ifftshift(shiftedSpectrum), cfg.nFFT) * sqrt(cfg.nFFT);
        symbolWithCp = [usefulSymbol(end-cfg.cpLen+1:end); usefulSymbol];
        outputIndices = (symbolIndex-1) * symbolLength + (1:symbolLength);
        waveform(outputIndices) = symbolWithCp;
    end
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
        case 'feedback'
            frameParams.infoBits = cfg.feedbackInfoBits;
            frameParams.codeRate = cfg.feedbackCodeRate;
            frameParams.scramblerSeed = cfg.scramblerSeedFeedback;
    end

    codedBitsPerOfdmSymbol = cfg.nDataSubcarriers * cfg.bitsPerSymbol;
    minCodedBits = ceil(frameParams.infoBits / frameParams.codeRate);
    frameParams.nDataSymbols = ceil(minCodedBits / codedBitsPerOfdmSymbol);
    frameParams.encodedBitCount = frameParams.nDataSymbols * codedBitsPerOfdmSymbol;
end
