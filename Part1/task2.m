clear; close all; clc;
% Parameters
Fs = 44100; % Sampling frequency
fc = 4000; % Carrier frequency
Tsym = 2.2676e-3; % Symbol time
Ns = round(Tsym * Fs); % Samples per symbol
beta = 2 + 2i; % Pilot
data = load('C:/Users/Hogwa/Desktop/SWE/4 PWC/Project_Part1/Signals_task2/Signal1.mat');
r = data.R(:);
N = length(r);
t = (0:N-1).' / Fs;
% Down Conversion
baseband = 2 * r .* exp(-1i * 2*pi*fc*t);
% Matched Filter
n = (0:Ns-1).';
pulse = sin(pi * n / (Ns-1));
pulse = pulse / sqrt(sum(abs(pulse).^2));
mf = flipud(pulse);
z = conv(baseband, mf, 'same');

figure
plot(abs(z));

% Synchronization
bestMetric =-inf;
bestPhase = 1;
for phase = 1:Ns
    samples = z(phase:Ns:end);
    metric = mean(abs(samples).^2);
    if metric > bestMetric
        bestMetric = metric;
        bestPhase = phase;
    end
end
symStream = z(bestPhase:Ns:end);
% Find Pilots
energy = abs(symStream).^2;
smoothEnergy = movmean(energy, 5);
threshold = 0.15 * max(smoothEnergy);
activeIdx = find(smoothEnergy > threshold);
startIdx = activeIdx(1);
endIdx = activeIdx(end);
rxSymbolsWithPilots = symStream(startIdx:endIdx);

figure
plot(abs(rxSymbolsWithPilots));

% Improve Pilot Detection
L = length(rxSymbolsWithPilots);
magSym = abs(rxSymbolsWithPilots);
pilotThreshold = 0.7 * max(magSym);
pilotCandidates = find(magSym > pilotThreshold);

if length(pilotCandidates) >= 2
    firstPilot = pilotCandidates(1);
    lastPilot = pilotCandidates(end);
    rxSymbolsWithPilots = rxSymbolsWithPilots(firstPilot:lastPilot);
end

% Channel Estimation
rxPilotStart = rxSymbolsWithPilots(1);
rxPilotEnd = rxSymbolsWithPilots(end);
h1 = rxPilotStart / beta;
h2 = rxPilotEnd / beta;
hHat = (h1 + h2) / 2;
% Pilot Removal
rxDataSymbols = rxSymbolsWithPilots(2:end-1);
% Channel Equalization
eqSymbols = rxDataSymbols / hHat;
% QPSK Decision
bits = [];
for k = 1:length(eqSymbols)
x = real(eqSymbols(k));
y = imag(eqSymbols(k));
if x >= 0 && y >= 0
bits = [bits 0 0];
elseif x < 0 && y >= 0
bits = [bits 1 0];
elseif x < 0 && y < 0
bits = [bits 1 1];
elseif x >= 0 && y < 0
bits = [bits 0 1];
end
end
% Decoding
TransmittedBits = bits;
wh = 2.^(6:-1:0);
m = char(TransmittedBits(1:7)*wh');
for l = 2:floor(length(TransmittedBits)/7)
m = [m char(TransmittedBits(7*(l-1)+1:7*l)*wh')];
end
fprintf('\nDecoded message:\n');
disp(m);

% plot
figure;
plot(real(z));
grid on;
title('After Downconversion and Matched Filter');
figure;
plot(abs(symStream));
grid on;
title('Right Phase');