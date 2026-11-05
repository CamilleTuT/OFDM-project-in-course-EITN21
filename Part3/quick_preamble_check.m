function [hasPacket, metricInfo] = quick_preamble_check(rxSignal, cfg, metricThreshold)
%QUICK_PREAMBLE_CHECK Fast preamble scan for Task 2 listening loops.

if nargin < 3
    metricThreshold = 0.20;
end

rxSignal = rxSignal(:);
[~, preambleInfo] = generate_preamble(cfg);

nFFT = cfg.nFFT;
cpLen = cfg.cpLen;
maxStart = numel(rxSignal) - cpLen - 2 * nFFT + 1;
if maxStart < 1
    hasPacket = false;
    metricInfo.peakMetric = 0;
    metricInfo.packetStartEstimate = NaN;
    metricInfo.reason = 'buffer too short';
    return;
end

pairProducts = conj(rxSignal(1:end-nFFT)) .* rxSignal(1+nFFT:end);
pairSums = filter(ones(nFFT, 1), 1, pairProducts);
pairSums = pairSums(nFFT:end);

sampleEnergy = abs(rxSignal).^2;
windowEnergy = filter(ones(nFFT, 1), 1, sampleEnergy);
windowEnergy = windowEnergy(nFFT:end);

k = (1:maxStart).';
pilotStart = k + cpLen;
firstEnergy = windowEnergy(pilotStart);
secondEnergy = windowEnergy(pilotStart + nFFT);
repeatedMetric = abs(pairSums(pilotStart)) ./ ...
    (sqrt(firstEnergy .* secondEnergy) + eps);

energyMetric = 0.5 * (firstEnergy + secondEnergy);
if max(energyMetric) > 0
    repeatedMetric(energyMetric < cfg.syncEnergyThreshold * max(energyMetric)) = 0;
end

knownCorr = conv(rxSignal, conj(flipud(preambleInfo.timePilot)), 'valid');
knownMetric = abs(knownCorr(pilotStart)) ./ ...
    (norm(preambleInfo.timePilot) * sqrt(firstEnergy) + eps);

combinedMetric = repeatedMetric .* knownMetric;
[peakMetric, peakIndex] = max(combinedMetric);

hasPacket = peakMetric >= metricThreshold;
metricInfo.peakMetric = peakMetric;
metricInfo.packetStartEstimate = peakIndex;
metricInfo.metricThreshold = metricThreshold;
metricInfo.repeatedMetricAtPeak = repeatedMetric(peakIndex);
metricInfo.knownMetricAtPeak = knownMetric(peakIndex);
metricInfo.reason = 'ok';
end
