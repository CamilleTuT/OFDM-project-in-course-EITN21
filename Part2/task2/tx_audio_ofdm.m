clear;
close all;

cfg = audio_config();
[packets, fileSummary] = make_packets(cfg);

fprintf('TX\n');
fprintf('File size: %d bytes, %d bits, %d packets.\n', ...
    fileSummary.numBytes, fileSummary.numBits, fileSummary.totalPackets);
fprintf('StartTime: %s\n', char(string(cfg.startTime, 'yyyy-MM-dd HH:mm:ss')));
fprintf('Round: data %.2fs, guard %.2fs, decode wait %.2fs, feedback %.2fs.\n', ...
    cfg.dataSlotSeconds, cfg.guardSeconds, cfg.decodeWaitSeconds, cfg.feedbackSlotSeconds);

pendingPacketIds = 1:fileSummary.totalPackets;
roundPeriodSeconds = cfg.dataSlotSeconds + cfg.guardSeconds + cfg.decodeWaitSeconds + cfg.feedbackSlotSeconds + cfg.guardSeconds;

for roundId = 1:cfg.maxRounds
    dataStartTime = cfg.startTime + seconds((roundId - 1) * roundPeriodSeconds);
    feedbackStartTime = dataStartTime + seconds(cfg.dataSlotSeconds + cfg.guardSeconds + cfg.decodeWaitSeconds);

    fprintf('\nTX round %d, pending packets: %d\n', roundId, numel(pendingPacketIds));

    dataAudio = build_data_round_audio(packets, pendingPacketIds, cfg);
    fprintf('Round %d sends %d packets: ', roundId, numel(pendingPacketIds));
    fprintf('%d ', pendingPacketIds);
    fprintf('\n');
    plot_audio_signal(dataAudio, cfg, sprintf('TX data round %d', roundId));

    wait_until_time(dataStartTime, sprintf('round %d data slot', roundId));
    player = audioplayer(dataAudio, cfg.audioFs, cfg.audioBitsPerSample);
    playblocking(player);

    wait_until_time(feedbackStartTime, sprintf('round %d feedback slot', roundId));
    fprintf('Recording feedback for %.2f s...\n', cfg.feedbackSlotSeconds);
    recorder = audiorecorder(cfg.audioFs, cfg.audioBitsPerSample, cfg.recordingChannels);
    recordblocking(recorder, cfg.feedbackSlotSeconds);
    feedbackAudio = getaudiodata(recorder, 'double');
    feedbackAudio = feedbackAudio(:);

    plot_audio_signal(feedbackAudio, cfg, sprintf('TX recorded feedback round %d', roundId));

    decodeTimer = tic;
    [feedbackFrames, ~] = decode_audio_stream(feedbackAudio, cfg, 'feedback');
    [feedback, validFeedbackCount] = parse_feedback_frames(feedbackFrames, cfg, roundId);
    decodeSeconds = toc(decodeTimer);

    fprintf('Feedback decode time: %.3f s, candidates: %d, valid decoded: %d.\n', ...
        decodeSeconds, numel(feedbackFrames), validFeedbackCount);

    if ~feedback.valid
        fprintf('No valid feedback decoded. Repeating the same pending set.\n');
        continue;
    end

    if feedback.status == cfg.feedbackStatusEndOfFile
        fprintf('Receiver reports END_OF_FILE. Transmission complete.\n');
        break;
    elseif feedback.status == cfg.feedbackStatusMissingBitmap
        pendingPacketIds = feedback.missingPacketIds;
        fprintf('Receiver requests %d packets for retransmission.\n', numel(pendingPacketIds));
    elseif feedback.status == cfg.feedbackStatusNoLock
        fprintf('Receiver reports NO_LOCK. Repeating the same pending set.\n');
    end
end

function [packets, fileSummary] = make_packets(cfg)
    fileId = fopen(cfg.inputFile, 'rb');
    bytes = fread(fileId, Inf, 'uint8=>uint8');
    fclose(fileId);

    fileBits = zeros(numel(bytes) * 8, 1);
    for k = 1:8
        fileBits(k:8:end) = double(bitget(bytes, 9-k));
    end

    totalPackets = ceil(numel(fileBits) / cfg.packetPayloadBits);

    packets = repmat(struct('infoBits', []), totalPackets, 1);

    for packetId = 1:totalPackets
        firstBit = (packetId - 1) * cfg.packetPayloadBits + 1;
        lastBit = min(packetId * cfg.packetPayloadBits, numel(fileBits));
        payloadBits = fileBits(firstBit:lastBit);
        payloadBitLength = numel(payloadBits);
        paddedPayload = [payloadBits; zeros(cfg.packetPayloadBits-payloadBitLength, 1)];

        packets(packetId).infoBits = build_data_packet_bits(cfg, ...
            cfg.packetTypeData, totalPackets, packetId, ...
            payloadBitLength, numel(fileBits), paddedPayload);
    end

    fileSummary.numBytes = numel(bytes);
    fileSummary.numBits = numel(fileBits);
    fileSummary.totalPackets = totalPackets;
end

function audio = build_data_round_audio(packets, pendingPacketIds, cfg)
    slotBasebandSamples = floor(cfg.dataSlotSeconds * cfg.basebandFs);
    gap = zeros(round(cfg.packetGapSeconds * cfg.basebandFs), 1);
    roundBaseband = [];

    for k = 1:numel(pendingPacketIds)
        packetId = pendingPacketIds(k);
        frame = build_phy_frame(packets(packetId).infoBits, cfg, 'data');
        roundBaseband = [roundBaseband; frame; gap];
    end

    roundBaseband = [roundBaseband; zeros(max(0, slotBasebandSamples-numel(roundBaseband)), 1)];
    roundBaseband = roundBaseband(1:slotBasebandSamples);

    audioComplex = resample(roundBaseband, cfg.audioFs, cfg.basebandFs);
    timeAxis = (0:numel(audioComplex)-1).' / cfg.audioFs;
    audio = real(audioComplex .* exp(1j * 2*pi*cfg.audioCarrierHz*timeAxis));
    peakAmplitude = max(abs(audio));
    audio = cfg.txAmplitude * audio / peakAmplitude;

    slotAudioSamples = floor(cfg.dataSlotSeconds * cfg.audioFs);
    audio = [audio(:); zeros(max(0, slotAudioSamples-numel(audio)), 1)];
    audio = audio(1:slotAudioSamples);
end

function bits = build_data_packet_bits(cfg, packetType, totalPackets, packetId, ...
        payloadBitLength, fileBitLength, payloadBits)
    payloadBits = double(payloadBits(:) ~= 0);
    payloadBits = [payloadBits; zeros(cfg.packetPayloadBits-numel(payloadBits), 1)];

    bits = [
        uint_to_bits(packetType, 4)
        uint_to_bits(totalPackets, 16)
        uint_to_bits(packetId, 16)
        uint_to_bits(payloadBitLength, 16)
        uint_to_bits(fileBitLength, 32)
        payloadBits
    ];
end

function [feedback, validFeedbackCount] = parse_feedback_frames(frames, cfg, expectedRoundId)
    feedback.valid = false;
    feedback.status = 0;
    feedback.missingPacketIds = [];
    validFeedbackCount = 0;

    for k = 1:numel(frames)
        if ~frames(k).valid
            continue;
        end
        validFeedbackCount = validFeedbackCount + 1;

        bits = double(frames(k).infoBits(:) ~= 0);
        bitPosition = 1;
        [status, bitPosition] = read_uint_bits(bits, bitPosition, 4);
        [roundId, bitPosition] = read_uint_bits(bits, bitPosition, 16);
        [totalPackets, bitPosition] = read_uint_bits(bits, bitPosition, 16);

        bitmap = logical(bits(bitPosition:bitPosition+cfg.maxFeedbackPackets-1));
        if double(roundId) == expectedRoundId
            feedback.valid = true;
            feedback.status = uint8(status);
            feedback.missingPacketIds = ...
                find(bitmap(1:double(totalPackets))).';
            return;
        end
    end
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
    remaining = seconds(targetTime - datetime('now'));
    fprintf('Waiting %.2f s until %s.\n', remaining, label);
    pause(remaining);
end

function plot_audio_signal(audio, cfg, titleText)
    audio = double(audio(:));

    timeAxis = (0:numel(audio)-1).' / cfg.audioFs;
    plotSamples = min(numel(audio), round(0.50 * cfg.audioFs));

    figure
    plot(timeAxis(1:plotSamples), audio(1:plotSamples));
    grid on;
    xlabel('Time (s)', 'FontSize', 18);
    ylabel('Amplitude', 'FontSize', 18);
    title([titleText ' waveform'], 'FontSize', 18);
    set(gca, 'FontSize', 18);

    spectrumFftLength = 2^nextpow2(max(4096, numel(audio)));
    window = 0.5 - 0.5*cos(2*pi*(0:numel(audio)-1).'/(numel(audio)-1));
    spectrum = fftshift(fft(audio .* window, spectrumFftLength));
    frequency = cfg.audioFs * (-spectrumFftLength/2:spectrumFftLength/2-1).' / spectrumFftLength;
    spectrumPeak = max(abs(spectrum));
    figure
    plot(frequency, 20*log10(abs(spectrum)/spectrumPeak + eps));
    grid on;
    xlabel('Frequency (Hz)', 'FontSize', 18);
    ylabel('Magnitude (dB, normalized)', 'FontSize', 18);
    title([titleText ' spectrum'], 'FontSize', 18);
    xlim([0 cfg.audioFs/2]);
    ylim([-120 5]);
    set(gca, 'FontSize', 18);
end
