function [hActive, debug] = estimate_channel_from_preamble(rxSignal, packetStart, cfg, preambleInfo)
%ESTIMATE_CHANNEL_FROM_PREAMBLE Initial frequency-domain channel estimate.

pilot1Index = packetStart + cfg.cpLen : packetStart + cfg.cpLen + cfg.nFFT - 1;
pilot2Index = packetStart + cfg.cpLen + cfg.nFFT : packetStart + cfg.cpLen + 2*cfg.nFFT - 1;
if pilot2Index(end) > numel(rxSignal)
    error('Received signal ended inside the preamble.');
end

pilotTime = 0.5 * (rxSignal(pilot1Index) + rxSignal(pilot2Index));
pilotFreqShifted = fftshift(fft(pilotTime, cfg.nFFT) / sqrt(cfg.nFFT));
rxActive = pilotFreqShifted(cfg.activeShiftIndices);

hActive = ones(cfg.nActive, 1);
nonZeroPilot = abs(preambleInfo.activePilot) > 0;
hActive(nonZeroPilot) = rxActive(nonZeroPilot) ./ preambleInfo.activePilot(nonZeroPilot);
hActive(~nonZeroPilot) = 1;

guardMask = true(cfg.nFFT, 1);
guardMask(cfg.activeShiftIndices) = false;
debug.rxActivePreamble = rxActive;
debug.hActive = hActive;
debug.guardNoisePower = mean(abs(pilotFreqShifted(guardMask)).^2);
end
