function cfg = audio_config()
    root = fileparts(mfilename('fullpath'));

    cfg.inputFile = fullfile(root, 'data', 'test1.txt');
    cfg.outputFile = fullfile(root, 'received', 'test1_received.txt');

    cfg.audioFs = 48000;
    cfg.basebandFs = 12000;
    cfg.audioCarrierHz = 10000;
    cfg.txAmplitude = 0.65;
    cfg.audioBitsPerSample = 16;
    cfg.recordingChannels = 1;

    cfg.nFFT = 128;
    cfg.nActive = 96;
    cfg.cpLen = 32;
    firstActiveFftIndex = cfg.nFFT/2 + 1 - cfg.nActive/2;
    cfg.activeFftIndices = (firstActiveFftIndex : firstActiveFftIndex + cfg.nActive - 1).';

    cfg.continuousPilotIndices = [12; 36; 61; 85];
    cfg.continuousPilotSymbols = [1; -1; 1; -1];
    cfg.dataActiveIndices = setdiff((1:cfg.nActive).', cfg.continuousPilotIndices, 'stable');
    cfg.nDataSubcarriers = numel(cfg.dataActiveIndices);

    cfg.modulation = 'QPSK';
    cfg.bitsPerSymbol = 2;
    cfg.dataCodeRate = 1/2;
    cfg.feedbackCodeRate = 1/3;
    cfg.ldpcMaxIterationsData = 25;
    cfg.ldpcMaxIterationsFeedback = 40;

    cfg.packetPayloadBits = 800;
    cfg.dataHeaderBits = 84;
    cfg.dataInfoBits = cfg.dataHeaderBits + cfg.packetPayloadBits;
    cfg.maxFeedbackPackets = 512;
    cfg.feedbackInfoBits = 4 + 16 + 16 + cfg.maxFeedbackPackets;

    cfg.packetTypeData = uint8(1);
    cfg.feedbackStatusMissingBitmap = uint8(1);
    cfg.feedbackStatusEndOfFile = uint8(2);
    cfg.feedbackStatusNoLock = uint8(3);

    cfg.scramblerSeedData = 100;
    cfg.scramblerSeedFeedback = 200;
    cfg.preambleSeed = 300;
    cfg.preambleAmplitude = 1;

    cfg.minimumChannelMagnitude = 1e-6;
    cfg.minNoiseVariance = 2e-4;

    cfg.startTime = datetime(2026, 6, 13, 17, 10, 30);
    cfg.dataSlotSeconds = 12.0;
    cfg.guardSeconds = 0.25;
    cfg.decodeWaitSeconds = 5.0;
    cfg.feedbackSlotSeconds = 2.0;
    cfg.maxRounds = 20;
    cfg.packetGapSeconds = 0.050;
    cfg.feedbackRepeatCount = 3;
    cfg.feedbackRepeatGapSeconds = 0.040;
end
