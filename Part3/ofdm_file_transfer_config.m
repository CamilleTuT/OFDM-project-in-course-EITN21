function cfg = ofdm_file_transfer_config()
%OFDM_FILE_TRANSFER_CONFIG Central parameter set for the Pluto OFDM link.

root = fileparts(mfilename('fullpath'));

cfg.inputFile = fullfile(root, 'data', 'test1.txt');
cfg.outputFile = fullfile(root, 'output', 'test1_recovered.txt');

cfg.centerFrequency = 2.45e9;
cfg.sampleRate = 20e6;
cfg.radioID = 'usb:0';
cfg.txGain = 0;
cfg.rxGain = 40;
cfg.txAmplitude = 0.55;
cfg.rxCaptureExtraSamples = 50000;
cfg.rxSettleTime = 0.20;

cfg.nFFT = 512;
cfg.nActive = 301;
cfg.cpLen = 128;
cfg.activeShiftIndices = (cfg.nFFT/2 + 1 - floor(cfg.nActive/2)) : ...
                         (cfg.nFFT/2 + 1 + floor(cfg.nActive/2));

cfg.contPilotActiveIndices = [1 101 201 301].';
cfg.dcActiveIndex = 151;
cfg.contPilotSymbols = ones(numel(cfg.contPilotActiveIndices), 1);
cfg.dataActiveIndices = setdiff((1:cfg.nActive).', ...
    [cfg.contPilotActiveIndices; cfg.dcActiveIndex], 'stable');
cfg.nDataSubcarriers = numel(cfg.dataActiveIndices);

cfg.payloadCodeRate = 1/2;
cfg.overheadCodeRate = 1/2;
cfg.overheadOfdmSymbols = 1;
cfg.ldpcMaxIterations = 50;

cfg.metadataMagic = uint32(hex2dec('0F0D'));
cfg.metadataVersion = uint8(1);
cfg.metadataBitLength = 140;

cfg.scramblerSeedOverhead = uint32(511);
cfg.scramblerSeedPayload = uint32(2027);
cfg.preambleSeed = uint32(1337);
cfg.preamblePilotAmplitude = 1;

cfg.syncThreshold = 0.75;
cfg.syncEnergyThreshold = 0.05;
cfg.maxFineTimingSearch = cfg.cpLen;
cfg.minNoiseVariance = 1e-5;
cfg.channelFloor = 1e-7;

cfg.enablePlots = true;
cfg.maxConstellationPoints = 12000;
end
