function [rxBits, rxInfo, debug] = decode_ofdm_file_frame(rxSignal, cfg)
%DECODE_OFDM_FILE_FRAME Detect, equalize, decode, and write one OFDM file packet.

rxSignal = rxSignal(:);
[~, preambleInfo] = generate_preamble(cfg);
[syncInfo, syncDebug] = detect_preamble_and_cfo(rxSignal, cfg, preambleInfo);
n = (0:numel(rxSignal)-1).';
rxCfoCorrected = rxSignal .* exp(-1j * 2 * pi * syncInfo.cfoNormalized * n);
[hActive, channelDebug] = estimate_channel_from_preamble( ...
    rxCfoCorrected, syncInfo.packetStart, cfg, preambleInfo);

symbolLength = cfg.nFFT + cfg.cpLen;
availableSymbols = floor((numel(rxCfoCorrected) - syncInfo.dataStart + 1) / symbolLength);
if availableSymbols < cfg.overheadOfdmSymbols
    error('Not enough received samples to decode the overhead.');
end

[rawActiveGrid, fullShiftedGrid] = waveform_to_active_grid( ...
    rxCfoCorrected, syncInfo.dataStart, availableSymbols, cfg);
[equalizedActiveGrid, pilotTrack, noiseVariance, equalizeDebug] = ...
    equalize_ofdm_symbols(rawActiveGrid, hActive, cfg);

modulation = 'QPSK';
bitsPerSymbol = 2;
codedBitsPerOfdmSymbol = cfg.nDataSubcarriers * bitsPerSymbol;
overheadEncodedBitCount = cfg.overheadOfdmSymbols * codedBitsPerOfdmSymbol;

overheadSymbols = reshape(equalizedActiveGrid(cfg.dataActiveIndices, ...
    1:cfg.overheadOfdmSymbols), [], 1);
overheadSymbols = overheadSymbols(1:overheadEncodedBitCount / bitsPerSymbol);
overheadLlrs = nrSymbolDemodulate(overheadSymbols, modulation, ...
    noiseVariance, 'DecisionType', 'soft');
[overheadDecodedScrambled, overheadCrcError] = ldpc_decode_bits( ...
    overheadLlrs, cfg.metadataBitLength, cfg.overheadCodeRate, ...
    modulation, cfg.ldpcMaxIterations);

if overheadCrcError ~= 0
    error('Overhead LDPC/CRC decode failed. Cannot know payload length safely.');
end

overheadDecodedBits = scramble_bits(overheadDecodedScrambled, cfg.scramblerSeedOverhead);
metadata = unpack_metadata(overheadDecodedBits, cfg);

requiredSymbols = cfg.overheadOfdmSymbols + metadata.payloadOfdmSymbols;
if availableSymbols < requiredSymbols
    error('Need %d data OFDM symbols, but only %d were captured.', ...
        requiredSymbols, availableSymbols);
end

payloadColumns = cfg.overheadOfdmSymbols + (1:metadata.payloadOfdmSymbols);
payloadSymbols = reshape(equalizedActiveGrid(cfg.dataActiveIndices, payloadColumns), [], 1);
payloadSymbols = payloadSymbols(1:metadata.payloadEncodedBitCount / bitsPerSymbol);
payloadLlrs = nrSymbolDemodulate(payloadSymbols, modulation, ...
    noiseVariance, 'DecisionType', 'soft');

[payloadDecodedScrambled, payloadCrcError] = ldpc_decode_bits(payloadLlrs, ...
    metadata.fileBitLength, cfg.payloadCodeRate, modulation, ...
    cfg.ldpcMaxIterations);
rxBits = scramble_bits(payloadDecodedScrambled, cfg.scramblerSeedPayload);
rxBits = rxBits(1:metadata.fileBitLength);

outputWritten = false;
if payloadCrcError == 0
    bits_to_file(rxBits, cfg.outputFile, metadata.fileByteLength);
    outputWritten = true;
end

rxInfo.syncInfo = syncInfo;
rxInfo.metadata = metadata;
rxInfo.noiseVariance = noiseVariance;
rxInfo.overheadCrcError = overheadCrcError;
rxInfo.payloadCrcError = payloadCrcError;
rxInfo.outputFile = cfg.outputFile;
rxInfo.outputWritten = outputWritten;
rxInfo.availableSymbols = availableSymbols;
rxInfo.requiredSymbols = requiredSymbols;

debug.rxSignal = rxSignal;
debug.rxCfoCorrected = rxCfoCorrected;
debug.sync = syncDebug;
debug.channel = channelDebug;
debug.fullShiftedGrid = fullShiftedGrid;
debug.rawActiveGrid = rawActiveGrid;
debug.equalizedActiveGrid = equalizedActiveGrid;
debug.equalize = equalizeDebug;
debug.pilotTrack = pilotTrack;
debug.rawDataSymbols = rawActiveGrid(cfg.dataActiveIndices, 1:requiredSymbols);
debug.equalizedDataSymbols = equalizedActiveGrid(cfg.dataActiveIndices, 1:requiredSymbols);
debug.rxInfo = rxInfo;

if cfg.enablePlots
    plot_ofdm_debug(debug, cfg, 'rx');
end
end
