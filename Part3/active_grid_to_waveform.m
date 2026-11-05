function waveform = active_grid_to_waveform(activeGrid, cfg)
%ACTIVE_GRID_TO_WAVEFORM OFDM modulate active sub-channel grids.

nSymbols = size(activeGrid, 2);
symbolLength = cfg.nFFT + cfg.cpLen;
waveform = zeros(nSymbols * symbolLength, 1);

for s = 1:nSymbols
    freqShifted = zeros(cfg.nFFT, 1);
    freqShifted(cfg.activeShiftIndices) = activeGrid(:, s);
    timeNoCp = ifft(ifftshift(freqShifted), cfg.nFFT) * sqrt(cfg.nFFT);
    timeWithCp = [timeNoCp(end-cfg.cpLen+1:end); timeNoCp];

    outIndex = (s-1) * symbolLength + (1:symbolLength);
    waveform(outIndex) = timeWithCp;
end
end
