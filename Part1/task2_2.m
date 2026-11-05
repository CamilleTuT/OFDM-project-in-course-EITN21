clear; clc; close all;

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
z = conv(baseband, mf);

% Synchronization
absz = abs(z);
[~, pos1] = max(absz);

tempdata = absz;
tempdata(pos1-Ns:pos1+Ns) = 0;
[~, pos2] = max(tempdata);

pilotPos = sort([pos1 pos2]);
startSample = pilotPos(1);
endSample = pilotPos(2);

numSym = round((endSample - startSample)/Ns) + 1;
sampleIdx = startSample + (0:numSym-1)*Ns;

rxSymbolsWithPilots = z(sampleIdx);

% Channel Estimation
rxPilotStart = rxSymbolsWithPilots(1);
rxPilotEnd   = rxSymbolsWithPilots(end);

h1 = rxPilotStart / beta;
h2 = rxPilotEnd   / beta;
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
title('After Downconversion and Matched Filter', 'FontSize', 18);
set(gca, 'FontSize', 18)
