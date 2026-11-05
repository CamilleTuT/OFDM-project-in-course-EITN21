function [equalizedActiveGrid, pilotTrack, noiseVariance, debug] = equalize_ofdm_symbols(rawActiveGrid, hActive, cfg)
%EQUALIZE_OFDM_SYMBOLS One-tap equalization plus continuous pilot CPE tracking.

hSafe = hActive(:);
hSafe(abs(hSafe) < cfg.channelFloor) = 1;
coarseEqualized = rawActiveGrid ./ repmat(hSafe, 1, size(rawActiveGrid, 2));

[pilotTrack, equalizedActiveGrid] = track_continuous_pilots(coarseEqualized, cfg);

pilotAfter = equalizedActiveGrid(cfg.contPilotActiveIndices, :);
expected = repmat(cfg.contPilotSymbols, 1, size(pilotAfter, 2));
pilotError = pilotAfter - expected;
noiseVariance = max(mean(abs(pilotError(:)).^2), cfg.minNoiseVariance);

debug.coarseEqualizedActiveGrid = coarseEqualized;
debug.pilotError = pilotError;
debug.noiseVariance = noiseVariance;
end
