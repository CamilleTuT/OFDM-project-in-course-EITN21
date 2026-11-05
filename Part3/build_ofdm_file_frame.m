function [txWaveform, txInfo, debug] = build_ofdm_file_frame(cfg)
%BUILD_OFDM_FILE_FRAME Convert one file into a complete OFDM packet waveform.

[fileBits, fileInfo] = file_to_bits(cfg.inputFile);
modulation = 'QPSK';
bitsPerSymbol = 2;
codedBitsPerOfdmSymbol = cfg.nDataSubcarriers * bitsPerSymbol;

payloadOfdmSymbols = ceil(ceil(numel(fileBits) / cfg.payloadCodeRate) / ...
    codedBitsPerOfdmSymbol);
payloadEncodedBitCount = payloadOfdmSymbols * codedBitsPerOfdmSymbol;

payloadInputBits = scramble_bits(fileBits, cfg.scramblerSeedPayload);
payloadCodedBits = ldpc_encode_bits(payloadInputBits, cfg.payloadCodeRate, ...
    modulation, payloadEncodedBitCount);
payloadSymbols = nrSymbolModulate(payloadCodedBits, modulation);

metadata.fileBitLength = numel(fileBits);
metadata.fileByteLength = fileInfo.numBytes;
metadata.payloadOfdmSymbols = payloadOfdmSymbols;
metadata.payloadEncodedBitCount = payloadEncodedBitCount;
overheadBits = pack_metadata(metadata, cfg);
overheadInputBits = scramble_bits(overheadBits, cfg.scramblerSeedOverhead);

overheadEncodedBitCount = cfg.overheadOfdmSymbols * codedBitsPerOfdmSymbol;
overheadCodedBits = ldpc_encode_bits(overheadInputBits, cfg.overheadCodeRate, ...
    modulation, overheadEncodedBitCount);
overheadSymbols = nrSymbolModulate(overheadCodedBits, modulation);

allDataSymbols = [overheadSymbols; payloadSymbols];
nOfdmDataSymbols = cfg.overheadOfdmSymbols + payloadOfdmSymbols;
expectedQamSymbols = nOfdmDataSymbols * cfg.nDataSubcarriers;
if numel(allDataSymbols) ~= expectedQamSymbols
    error('QAM symbol count mismatch while building the frame.');
end

activeGrid = zeros(cfg.nActive, nOfdmDataSymbols);
for s = 1:nOfdmDataSymbols
    active = zeros(cfg.nActive, 1);
    active(cfg.contPilotActiveIndices) = cfg.contPilotSymbols;
    active(cfg.dcActiveIndex) = 0;
    dataIndex = (s-1) * cfg.nDataSubcarriers + (1:cfg.nDataSubcarriers);
    active(cfg.dataActiveIndices) = allDataSymbols(dataIndex);
    activeGrid(:, s) = active;
end

[preamble, preambleInfo] = generate_preamble(cfg);
dataWaveform = active_grid_to_waveform(activeGrid, cfg);
txWaveform = [preamble; dataWaveform];

peakAmplitude = max(abs(txWaveform));
if peakAmplitude > 0
    txWaveform = cfg.txAmplitude * txWaveform / peakAmplitude;
end

txInfo.fileInfo = fileInfo;
txInfo.fileBitLength = numel(fileBits);
txInfo.fileByteLength = fileInfo.numBytes;
txInfo.payloadOfdmSymbols = payloadOfdmSymbols;
txInfo.payloadEncodedBitCount = payloadEncodedBitCount;
txInfo.overheadEncodedBitCount = overheadEncodedBitCount;
txInfo.totalOfdmDataSymbols = nOfdmDataSymbols;
txInfo.waveformSamples = numel(txWaveform);
txInfo.waveformDurationSeconds = numel(txWaveform) / cfg.sampleRate;

debug.txWaveform = txWaveform;
debug.activeGrid = activeGrid;
debug.preambleInfo = preambleInfo;
debug.txInfo = txInfo;

if cfg.enablePlots
    plot_ofdm_debug(debug, cfg, 'tx');
end
end
