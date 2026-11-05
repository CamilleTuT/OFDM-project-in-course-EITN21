function [pilotTrack, correctedActiveGrid] = track_continuous_pilots(equalizedActiveGrid, cfg)
%TRACK_CONTINUOUS_PILOTS Estimate and remove common residual phase per symbol.

pilotRx = equalizedActiveGrid(cfg.contPilotActiveIndices, :);
expected = repmat(cfg.contPilotSymbols, 1, size(pilotRx, 2));
pilotRatio = pilotRx ./ expected;
commonRatio = mean(pilotRatio, 1);

bad = abs(commonRatio) < cfg.channelFloor;
commonRatio(bad) = 1;

correctedActiveGrid = equalizedActiveGrid ./ repmat(commonRatio, cfg.nActive, 1);

pilotTrack.amplitude = abs(pilotRx);
pilotTrack.phase = unwrap(angle(pilotRatio), [], 2);
pilotTrack.commonAmplitude = abs(commonRatio);
pilotTrack.commonPhase = unwrap(angle(commonRatio));
pilotTrack.commonRatio = commonRatio;
end
