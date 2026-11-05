clear;
close all;

cfg = ofdm_file_transfer_config();
cfg.enablePlots = false;

leadTimeSeconds = 3;
repeatWindowSeconds = 60;
interPacketGapSeconds = 0.002;

[txWaveform, ~] = build_ofdm_file_frame(cfg);
interPacketGap = zeros(round(interPacketGapSeconds * cfg.sampleRate), 1);
repeatWaveform = [txWaveform; interPacketGap];

fprintf('TX starts after %.1f s and repeats for %.1f s.\n', ...
    leadTimeSeconds, repeatWindowSeconds);
pause(leadTimeSeconds);

tx = sdrtx('Pluto', 'RadioID', cfg.radioID);
tx.CenterFrequency = cfg.centerFrequency;
tx.BasebandSampleRate = cfg.sampleRate;
tx.Gain = cfg.txGain;

cleanupTx = onCleanup(@() release(tx));
transmitRepeat(tx, repeatWaveform);

fprintf('TX repeat is active now.\n');
pause(repeatWindowSeconds);
fprintf('TX repeat window complete.\n');
