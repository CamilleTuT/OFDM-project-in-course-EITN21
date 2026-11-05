clear;
close all;

cfg = ofdm_file_transfer_config();

[txWaveform, txInfo] = build_ofdm_file_frame(cfg);
fprintf('Prepared %s: %d bytes, %d samples, %.3f ms.\n', ...
    txInfo.fileInfo.name, txInfo.fileByteLength, txInfo.waveformSamples, ...
    1e3 * txInfo.waveformDurationSeconds);

tx = sdrtx('Pluto', 'RadioID', cfg.radioID);
tx.CenterFrequency = cfg.centerFrequency;
tx.BasebandSampleRate = cfg.sampleRate;
tx.Gain = cfg.txGain;

rx = sdrrx('Pluto', 'RadioID', cfg.radioID);
rx.CenterFrequency = cfg.centerFrequency;
rx.BasebandSampleRate = cfg.sampleRate;
rx.GainSource = 'Manual';
rx.Gain = cfg.rxGain;
rx.OutputDataType = 'double';
rx.SamplesPerFrame = numel(txWaveform) + cfg.rxCaptureExtraSamples;

transmitRepeat(tx, txWaveform);
pause(cfg.rxSettleTime);
rxSignal = rx();

release(tx);
release(rx);

[~, rxInfo] = decode_ofdm_file_frame(rxSignal, cfg);
fprintf('Decoded output: %s\n', rxInfo.outputFile);
fprintf('Payload CRC error flag: %d\n', rxInfo.payloadCrcError);
