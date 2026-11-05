function [activeGrid, fullShiftedGrid] = waveform_to_active_grid(rxSignal, dataStart, nSymbols, cfg)
%WAVEFORM_TO_ACTIVE_GRID Remove CP and FFT OFDM data symbols.

rxSignal = rxSignal(:);
symbolLength = cfg.nFFT + cfg.cpLen;
activeGrid = zeros(cfg.nActive, nSymbols);
fullShiftedGrid = zeros(cfg.nFFT, nSymbols);

for s = 1:nSymbols
    cpStart = dataStart + (s-1) * symbolLength;
    fftStart = cpStart + cfg.cpLen;
    fftIndex = fftStart : fftStart + cfg.nFFT - 1;
    if fftIndex(end) > numel(rxSignal)
        error('Received signal ended inside OFDM symbol %d.', s);
    end

    timeNoCp = rxSignal(fftIndex);
    freqShifted = fftshift(fft(timeNoCp, cfg.nFFT) / sqrt(cfg.nFFT));
    fullShiftedGrid(:, s) = freqShifted;
    activeGrid(:, s) = freqShifted(cfg.activeShiftIndices);
end
end
