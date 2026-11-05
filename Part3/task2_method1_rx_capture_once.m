clear;
close all;

cfg = ofdm_file_transfer_config();
cfg.enablePlots = true;

captureSeconds = 0.050;
rxSamplesPerFrame = ceil(captureSeconds * cfg.sampleRate);

rx = sdrrx('Pluto', 'RadioID', cfg.radioID);
rx.CenterFrequency = cfg.centerFrequency;
rx.BasebandSampleRate = cfg.sampleRate;
rx.GainSource = 'Manual';
rx.Gain = cfg.rxGain;
rx.OutputDataType = 'double';
rx.SamplesPerFrame = rxSamplesPerFrame;

cleanupRx = onCleanup(@() release(rx));
rxSignal = rx();

[hasPacket, metricInfo] = quick_preamble_check(rxSignal, cfg, 0.20);
fprintf('Preamble metric: %.3f, candidate start sample: %g.\n', ...
    metricInfo.peakMetric, metricInfo.packetStartEstimate);
if ~hasPacket
    error('No strong preamble was found. Capture again during the TX repeat window.');
end

[~, rxInfo] = decode_ofdm_file_frame(rxSignal, cfg);
fprintf('Decoded output: %s\n', rxInfo.outputFile);
fprintf('Overhead CRC error flag: %d\n', rxInfo.overheadCrcError);
fprintf('Payload CRC error flag: %d\n', rxInfo.payloadCrcError);
fprintf('Estimated CFO: %.2f Hz\n', rxInfo.syncInfo.cfoHz);
