function [preamble, info] = generate_preamble(cfg)
%GENERATE_PREAMBLE Build [CP][pilot][pilot] with identical pilot symbols.

activePilot = zeros(cfg.nActive, 1);
nonDcActive = setdiff((1:cfg.nActive).', cfg.dcActiveIndex, 'stable');
pilotBits = prbs_bits(numel(nonDcActive), cfg.preambleSeed);
pilotSigns = 2 * pilotBits - 1;
activePilot(nonDcActive) = cfg.preamblePilotAmplitude * pilotSigns;
activePilot(cfg.dcActiveIndex) = 0;

freqShifted = zeros(cfg.nFFT, 1);
freqShifted(cfg.activeShiftIndices) = activePilot;
timePilot = ifft(ifftshift(freqShifted), cfg.nFFT) * sqrt(cfg.nFFT);

preamble = [timePilot(end-cfg.cpLen+1:end); timePilot; timePilot];

info.activePilot = activePilot;
info.freqShifted = freqShifted;
info.timePilot = timePilot;
info.length = numel(preamble);
end
