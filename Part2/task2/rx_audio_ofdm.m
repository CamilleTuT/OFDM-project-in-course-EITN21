clear;
close all;

cfg = audio_config();
[referencePackets, ~] = make_packets(cfg);

fprintf('RX\n');
packetStore = {};
totalPackets = 0;
fileBitLength = 0;
missingPacketIds = [];
roundPeriodSeconds = cfg.dataSlotSeconds + cfg.guardSeconds + cfg.decodeWaitSeconds + cfg.feedbackSlotSeconds + cfg.guardSeconds;

for roundId = 1:cfg.maxRounds
    dataStartTime = cfg.startTime + seconds((roundId - 1) * roundPeriodSeconds);
    feedbackStartTime = dataStartTime + seconds(cfg.dataSlotSeconds + cfg.guardSeconds + cfg.decodeWaitSeconds);

    fprintf('\nRX round %d\n', roundId);
    wait_until_time(dataStartTime, sprintf('round %d data slot', roundId));

    fprintf('Recording data for %.2f s...\n', cfg.dataSlotSeconds);
    recorder = audiorecorder(cfg.audioFs, cfg.audioBitsPerSample, cfg.recordingChannels);
    recordblocking(recorder, cfg.dataSlotSeconds);
    dataAudio = getaudiodata(recorder, 'double');
    dataAudio = dataAudio(:);

    decodeTimer = tic;
    [dataFrames, syncMetric] = decode_audio_stream(dataAudio, cfg, 'data');
    decodedPackets = parse_data_frames(dataFrames, cfg);
    decodeSeconds = toc(decodeTimer);
    fprintf('Data decode time: %.3f s, candidates: %d, valid packets: %d.\n', decodeSeconds, numel(dataFrames), numel(decodedPackets));

    analyze_rx_round(dataAudio, dataFrames, syncMetric, decodedPackets, referencePackets, cfg, roundId);

    for k = 1:numel(decodedPackets)
        packet = decodedPackets(k);
        if totalPackets == 0
            totalPackets = packet.totalPackets;
            fileBitLength = packet.fileBitLength;
            packetStore = cell(totalPackets, 1);
        end
        if packet.totalPackets == totalPackets && packet.packetId >= 1 && packet.packetId <= totalPackets
            packetStore{packet.packetId} = packet.payloadBits;
        end
    end

    if totalPackets == 0
        missingPacketIds = [];
        status = cfg.feedbackStatusNoLock;
        fprintf('No valid data packet yet; sending NO_LOCK feedback.\n');
    else
        missingPacketIds = find(cellfun(@isempty, packetStore)).';
        if isempty(missingPacketIds)
            status = cfg.feedbackStatusEndOfFile;
            fprintf('All %d packets decoded; sending END_OF_FILE feedback.\n', totalPackets);
        else
            status = cfg.feedbackStatusMissingBitmap;
            fprintf('%d/%d packets still missing.\n', numel(missingPacketIds), totalPackets);
        end
    end

    feedbackAudio = build_feedback_audio(missingPacketIds, roundId, totalPackets, cfg, status);
    plot_audio_signal(feedbackAudio, cfg, sprintf('RX feedback round %d', roundId));

    wait_until_time(feedbackStartTime, sprintf('round %d feedback slot', roundId));
    player = audioplayer(feedbackAudio, cfg.audioFs, cfg.audioBitsPerSample);
    playblocking(player);

    if totalPackets > 0 && isempty(missingPacketIds)
        reconstruct_file(packetStore, fileBitLength, cfg.outputFile);
        compare_recovered_file(cfg.inputFile, cfg.outputFile, roundId, cfg);
        break;
    end
end

function [packets, fileSummary] = make_packets(cfg)
    [fileBits, fileSummary] = file_to_bits(cfg.inputFile);
    totalPackets = ceil(numel(fileBits) / cfg.packetPayloadBits);

    packets = repmat(struct('infoBits', []), totalPackets, 1);

    for packetId = 1:totalPackets
        firstBit = (packetId - 1) * cfg.packetPayloadBits + 1;
        lastBit = min(packetId * cfg.packetPayloadBits, numel(fileBits));
        payloadBits = fileBits(firstBit:lastBit);
        payloadBitLength = numel(payloadBits);
        paddedPayload = [payloadBits; zeros(cfg.packetPayloadBits-payloadBitLength, 1)];

        packets(packetId).infoBits = [
            uint_to_bits(cfg.packetTypeData, 4)
            uint_to_bits(totalPackets, 16)
            uint_to_bits(packetId, 16)
            uint_to_bits(payloadBitLength, 16)
            uint_to_bits(numel(fileBits), 32)
            paddedPayload
        ];
    end

    fileSummary.totalPackets = totalPackets;
end

function decodedPackets = parse_data_frames(frames, cfg)
    decodedPackets = repmat(struct('packetId', 0, 'totalPackets', 0, ...
        'payloadBitLength', 0, 'fileBitLength', 0, ...
        'payloadBits', [], 'sourceFrameIndex', 0), 0, 1);

    for k = 1:numel(frames)
        if ~frames(k).valid
            continue;
        end

        bits = double(frames(k).infoBits(:) ~= 0);
        bitPosition = 1;
        [packetType, bitPosition] = read_uint_bits(bits, bitPosition, 4);
        [totalPackets, bitPosition] = read_uint_bits(bits, bitPosition, 16);
        [packetId, bitPosition] = read_uint_bits(bits, bitPosition, 16);
        [payloadBitLength, bitPosition] = read_uint_bits(bits, bitPosition, 16);
        [fileBitLength, bitPosition] = read_uint_bits(bits, bitPosition, 32);

        if uint8(packetType) == cfg.packetTypeData && packetId >= 1 && packetId <= totalPackets
            paddedPayload = bits(bitPosition:bitPosition+cfg.packetPayloadBits-1);
            decodedPackets(end+1).packetId = double(packetId);
            decodedPackets(end).totalPackets = double(totalPackets);
            decodedPackets(end).payloadBitLength = double(payloadBitLength);
            decodedPackets(end).fileBitLength = double(fileBitLength);
            decodedPackets(end).payloadBits = paddedPayload(1:double(payloadBitLength));
            decodedPackets(end).sourceFrameIndex = k;
        end
    end
end

function audio = build_feedback_audio(missingPacketIds, roundId, totalPackets, cfg, status)
    bitmap = false(cfg.maxFeedbackPackets, 1);
    bitmap(missingPacketIds) = true;

    feedbackBits = [
        uint_to_bits(status, 4)
        uint_to_bits(roundId, 16)
        uint_to_bits(totalPackets, 16)
        double(bitmap)
    ];

    feedbackFrame = build_phy_frame(feedbackBits, cfg, 'feedback');
    gap = zeros(round(cfg.feedbackRepeatGapSeconds*cfg.basebandFs), 1);

    baseband = [];
    for k = 1:cfg.feedbackRepeatCount
        baseband = [baseband; feedbackFrame; gap];
    end

    slotBasebandSamples = floor(cfg.feedbackSlotSeconds*cfg.basebandFs);
    baseband = [baseband; zeros(max(0, slotBasebandSamples-numel(baseband)), 1)];
    baseband = baseband(1:slotBasebandSamples);

    audioComplex = resample(baseband, cfg.audioFs, cfg.basebandFs);
    timeAxis = (0:numel(audioComplex)-1).' / cfg.audioFs;
    audio = real(audioComplex .* exp(1j*2*pi*cfg.audioCarrierHz*timeAxis));
    peakAmplitude = max(abs(audio));
    audio = cfg.txAmplitude * audio / peakAmplitude;

    slotAudioSamples = floor(cfg.feedbackSlotSeconds*cfg.audioFs);
    audio = [audio(:); zeros(max(0, slotAudioSamples-numel(audio)), 1)];
    audio = audio(1:slotAudioSamples);
end

function reconstruct_file(packetStore, fileBitLength, outputFile)
    fileBits = vertcat(packetStore{:});
    fileBits = fileBits(1:fileBitLength);
    numBytes = fileBitLength/8;
    bits = fileBits ~= 0;
    bytes = zeros(numBytes, 1, 'uint8');
    for k = 1:8
        bytes = bitor(bytes, uint8(bits(k:8:end))*uint8(2^(8-k)));
    end

    fileId = fopen(outputFile, 'wb');
    fwrite(fileId, bytes, 'uint8');
    fclose(fileId);
end

function compare_recovered_file(inputFile, outputFile, completedRounds, cfg)
    [inputBits, inputInfo] = file_to_bits(inputFile);
    [outputBits, outputInfo] = file_to_bits(outputFile);
    comparedBitCount = min(numel(inputBits), numel(outputBits));
    bitErrors = sum(inputBits(1:comparedBitCount) ~= outputBits(1:comparedBitCount)) + abs(numel(inputBits)-numel(outputBits));

    roundPeriodSeconds = cfg.dataSlotSeconds + cfg.guardSeconds + cfg.decodeWaitSeconds + cfg.feedbackSlotSeconds + cfg.guardSeconds;
    elapsedSeconds = completedRounds*roundPeriodSeconds;

    endToEndPayloadRateBps = numel(inputBits)/elapsedSeconds;
    dataSlotOnlyPayloadRateBps = numel(inputBits) / (completedRounds*cfg.dataSlotSeconds);

    fprintf('\nRecovered file comparison\n');
    fprintf('Original file: (%d bytes, %d bits)\n', inputInfo.numBytes, numel(inputBits));
    fprintf('Recovered file: (%d bytes, %d bits)\n', outputInfo.numBytes, numel(outputBits));
    fprintf('Final decoded bit errors after LDPC/ARQ: %d\n', bitErrors);
    fprintf('End-to-end practical payload rate: %.2f bit/s\n', endToEndPayloadRateBps);
    fprintf('Data-slot-only practical payload rate: %.2f bit/s\n', dataSlotOnlyPayloadRateBps);
end

function analyze_rx_round(dataAudio, dataFrames, syncMetric, decodedPackets, packets, cfg, roundId)
    plot_audio_signal(dataAudio, cfg, sprintf('RX recorded data round %d', roundId));

    figure
    syncTimeAxis = (0:numel(syncMetric)-1).' / cfg.basebandFs;
    plot(syncTimeAxis, syncMetric);
    grid on;
    xlabel('Baseband time / s', 'FontSize', 18);
    ylabel('Preamble correlation', 'FontSize', 18);
    title(sprintf('Preamble synchronization metric, round %d', roundId), 'FontSize', 18);
    set(gca, 'FontSize', 18);

    receivedDataSymbols = [];
    phaseBeforeCfo = [];
    phaseAfterCompensation = [];
    for k = 1:numel(dataFrames)
        if ~isempty(dataFrames(k).analysis)
            receivedDataSymbols = [receivedDataSymbols; dataFrames(k).analysis.equalizedDataSymbols(:)];
            phaseBeforeCfo = [phaseBeforeCfo, dataFrames(k).analysis.pilotPhaseBeforeCfo, nan(4, 1)];
            phaseAfterCompensation = [phaseAfterCompensation, ...
                dataFrames(k).analysis.pilotPhaseAfterCompensation, nan(4, 1)];
        end
    end

    figure
    plot(real(receivedDataSymbols), imag(receivedDataSymbols), '.');
    hold on;
    idealQpsk = nrSymbolModulate([0; 0; 0; 1; 1; 0; 1; 1], 'QPSK');
    plot(real(idealQpsk), imag(idealQpsk), 'rx');
    grid on;
    axis equal;
    xlabel('In-phase', 'FontSize', 18);
    ylabel('Quadrature', 'FontSize', 18);
    title('Received equalized QPSK constellation', 'FontSize', 18);
    legend('Received symbols', 'Ideal QPSK points', 'FontSize', 16);
    set(gca, 'FontSize', 18);

    figure
    plot(phaseBeforeCfo.');
    grid on;
    xlabel('OFDM symbol number', 'FontSize', 18);
    ylabel('Pilot phase / rad', 'FontSize', 18);
    title(sprintf('Pilot phase before CFO compensation, round %d', roundId), 'FontSize', 18);
    legend('Pilot 1', 'Pilot 2', 'Pilot 3', 'Pilot 4', 'FontSize', 16);
    set(gca, 'FontSize', 18);

    figure
    plot(phaseAfterCompensation.');
    grid on;
    xlabel('OFDM symbol number', 'FontSize', 18);
    ylabel('Pilot phase / rad', 'FontSize', 18);
    title(sprintf('Pilot phase after all compensation, round %d', roundId), 'FontSize', 18);
    legend('Pilot 1', 'Pilot 2', 'Pilot 3', 'Pilot 4', 'FontSize', 16);
    set(gca, 'FontSize', 18);

    bitErrors = 0;
    bitCount = 0;
    for k = 1:numel(decodedPackets)
        packet = decodedPackets(k);
        decodedFrame = dataFrames(packet.sourceFrameIndex);

        [~, txCodedBits] = build_phy_frame(packets(packet.packetId).infoBits, cfg, 'data');
        hardBits = nrSymbolDemodulate(decodedFrame.analysis.equalizedDataSymbols, 'QPSK', 'DecisionType', 'hard');
        txBits = double(txCodedBits(:) ~= 0);
        comparedBitCount = min(numel(txBits), numel(hardBits));
        bitErrors = bitErrors + sum(txBits(1:comparedBitCount) ~= hardBits(1:comparedBitCount));
        bitCount = bitCount + comparedBitCount;
    end

    if bitCount == 0
        fprintf('Raw BER before LDPC decoding unavailable for this round.\n');
    else
        rawBer = bitErrors/bitCount;
        fprintf('Raw BER before LDPC decoding, decoded data frames only: %.6g (%d/%d bits)\n', rawBer, bitErrors, bitCount);
    end

    symbolSeconds = (cfg.nFFT+cfg.cpLen)/cfg.basebandFs;
    dataInfoBits = cfg.dataHeaderBits + cfg.packetPayloadBits;
    minCodedBits = ceil(dataInfoBits/cfg.dataCodeRate);
    codedBitsPerSymbol = cfg.nDataSubcarriers*cfg.bitsPerSymbol;
    nDataSymbols = ceil(minCodedBits/codedBitsPerSymbol);
    frameSeconds = (cfg.cpLen + 2*cfg.nFFT + nDataSymbols*(cfg.nFFT+cfg.cpLen))/cfg.basebandFs;

    rawCodedPhyRate = cfg.nDataSubcarriers*cfg.bitsPerSymbol/symbolSeconds;
    theoreticalPayloadRate = cfg.packetPayloadBits/(frameSeconds+cfg.packetGapSeconds);
    practicalPayloadRate = sum([decodedPackets.payloadBitLength])/cfg.dataSlotSeconds;

    fprintf('Raw coded PHY rate on data subcarriers: %.2f bit/s\n', rawCodedPhyRate);
    fprintf('Theoretical payload rate per full packet: %.2f bit/s\n', theoreticalPayloadRate);
    fprintf('This-round decoded payload rate during data slot: %.2f bit/s\n', practicalPayloadRate);
end

function [bits, fileSummary] = file_to_bits(filePath)
    fileId = fopen(filePath, 'rb');
    bytes = fread(fileId, Inf, 'uint8=>uint8');
    fclose(fileId);

    bits = zeros(numel(bytes)*8, 1);
    for k = 1:8
        bits(k:8:end) = double(bitget(bytes, 9-k));
    end

    fileSummary.numBytes = numel(bytes);
end

function bits = uint_to_bits(value, width)
    value = uint64(value);
    bits = zeros(width, 1);
    for k = 1:width
        bits(k) = double(bitget(value, width-k+1));
    end
end

function [value, nextBitPosition] = read_uint_bits(bits, bitPosition, width)
    value = uint64(0);
    for k = bitPosition:bitPosition+width-1
        value = bitshift(value, 1) + uint64(bits(k));
    end
    nextBitPosition = bitPosition + width;
end

function wait_until_time(targetTime, label)
    remaining = seconds(targetTime-datetime('now'));
    fprintf('Waiting %.2f s until %s.\n', remaining, label);
    pause(remaining);
end

function plot_audio_signal(audio, cfg, titleText)
    audio = double(audio(:));

    timeAxis = (0:numel(audio)-1).' / cfg.audioFs;
    plotSamples = min(numel(audio), round(0.50*cfg.audioFs));
    figure
    plot(timeAxis(1:plotSamples), audio(1:plotSamples));
    grid on;
    xlabel('Time (s)', 'FontSize', 18);
    ylabel('Amplitude', 'FontSize', 18);
    title([titleText ' waveform'], 'FontSize', 18);
    set(gca, 'FontSize', 18);

    spectrumFftLength = 2^nextpow2(max(4096, numel(audio)));
    window = 0.5 - 0.5*cos( ...
        2*pi*(0:numel(audio)-1).'/(numel(audio)-1));
    spectrum = fftshift(fft(audio.*window, spectrumFftLength));
    frequency = cfg.audioFs*(-spectrumFftLength/2:spectrumFftLength/2-1).'/spectrumFftLength;
    spectrumPeak = max(abs(spectrum));
    figure
    plot(frequency, 20*log10(abs(spectrum)/spectrumPeak+eps));
    grid on;
    xlabel('Frequency (Hz)', 'FontSize', 18);
    ylabel('Magnitude (dB, normalized)', 'FontSize', 18);
    title([titleText ' spectrum'], 'FontSize', 18);
    xlim([0 cfg.audioFs/2]);
    ylim([-120 5]);
    set(gca, 'FontSize', 18);
end
