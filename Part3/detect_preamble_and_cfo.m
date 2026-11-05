function [syncInfo, debug] = detect_preamble_and_cfo(rxSignal, cfg, preambleInfo)
%DETECT_PREAMBLE_AND_CFO Packet detection and coarse CFO from repeated pilot.

rxSignal = rxSignal(:);
maxStart = numel(rxSignal) - cfg.cpLen - 2 * cfg.nFFT + 1;
if maxStart < 1
    error('Received signal is shorter than one preamble.');
end

metric = zeros(maxStart, 1);
energyMetric = zeros(maxStart, 1);

for k = 1:maxStart
    firstPilot = rxSignal(k+cfg.cpLen : k+cfg.cpLen+cfg.nFFT-1);
    secondPilot = rxSignal(k+cfg.cpLen+cfg.nFFT : k+cfg.cpLen+2*cfg.nFFT-1);
    corrValue = sum(conj(firstPilot) .* secondPilot);
    firstEnergy = sum(abs(firstPilot).^2);
    secondEnergy = sum(abs(secondPilot).^2);
    energy = sqrt(firstEnergy * secondEnergy) + eps;
    energyMetric(k) = 0.5 * (firstEnergy + secondEnergy);
    metric(k) = abs(corrValue) / energy;
end

energyMask = energyMetric >= cfg.syncEnergyThreshold * max(energyMetric);
metric(~energyMask) = 0;

knownMetricAll = zeros(maxStart, 1);
for k = 1:maxStart
    firstPilot = rxSignal(k+cfg.cpLen : k+cfg.cpLen+cfg.nFFT-1);
    knownMetricAll(k) = abs(sum(conj(preambleInfo.timePilot) .* firstPilot)) / ...
        (norm(preambleInfo.timePilot) * norm(firstPilot) + eps);
end

combinedMetric = metric .* knownMetricAll;
if max(combinedMetric) == 0
    combinedMetric = metric;
end

[peakMetric, peakIndex] = max(combinedMetric);
candidates = find(combinedMetric >= cfg.syncThreshold * peakMetric);
if isempty(candidates)
    coarseStart = peakIndex;
else
    coarseStart = candidates(1);
end

searchStart = max(1, coarseStart - cfg.maxFineTimingSearch);
searchEnd = min(maxStart, coarseStart + cfg.maxFineTimingSearch);
knownMetric = zeros(searchEnd - searchStart + 1, 1);
for k = searchStart:searchEnd
    firstPilot = rxSignal(k+cfg.cpLen : k+cfg.cpLen+cfg.nFFT-1);
    knownMetric(k - searchStart + 1) = abs(sum(conj(preambleInfo.timePilot) .* firstPilot)) / ...
        (norm(preambleInfo.timePilot) * norm(firstPilot) + eps);
end
[~, fineOffset] = max(knownMetric);
packetStart = searchStart + fineOffset - 1;

firstPilot = rxSignal(packetStart+cfg.cpLen : packetStart+cfg.cpLen+cfg.nFFT-1);
secondPilot = rxSignal(packetStart+cfg.cpLen+cfg.nFFT : packetStart+cfg.cpLen+2*cfg.nFFT-1);
corrAtStart = sum(conj(firstPilot) .* secondPilot);
cfoNormalized = angle(corrAtStart) / (2 * pi * cfg.nFFT);

syncInfo.packetStart = packetStart;
syncInfo.dataStart = packetStart + cfg.cpLen + 2 * cfg.nFFT;
syncInfo.cfoNormalized = cfoNormalized;
syncInfo.cfoHz = cfoNormalized * cfg.sampleRate;
syncInfo.peakMetric = peakMetric;
syncInfo.coarseStart = coarseStart;

debug.correlationMetric = combinedMetric;
debug.repeatedCorrelationMetric = metric;
debug.knownPreambleMetricAll = knownMetricAll;
debug.correlationEnergy = energyMetric;
debug.correlationPeakIndex = peakIndex;
debug.knownPreambleMetric = knownMetric;
debug.knownPreambleSearchStart = searchStart;
end
