clear; clc; close all;

Fs = 44100;
fc = 10000;

Nsc = 128;
Ncp = 20;
Tsym = 58e-3;

M = round(Fs * Tsym / Nsc);
Lsym = Nsc + Ncp;

signalFile = 'C:/Users/Hogwa/Desktop/SWE/4 PWC/Project_Part1/Signals_task3/signal8.mat';

data = load(signalFile);
r = data.R(:);

t = (0:length(r)-1).' / Fs;

% Down conversion
baseband = 2 * r .* exp(-1i * 2*pi*fc*t);

% Low-pass filter
[B, A] = butter(8, 0.05);
yFiltered = filter(B, A, baseband);

[H, f] = freqz(B, A, 2048, Fs);
figure;
plot(f, abs(H));
grid on;
xlabel('Frequency (Hz)', 'FontSize', 18);
ylabel('Magnitude', 'FontSize', 18);
title('Amplitude Frequency Response of the Low-pass Filter', 'FontSize', 18);
set(gca, 'FontSize', 18)

% Generate pilot
randn('state', 100);
P = sign(randn(1, Nsc/2));

pilotTx = zeros(Nsc, 1);
pilotTx(1:2:end) = 2 * P(:);

pilotIndex = 1:2:Nsc;

% Convolutional code
trellis = poly2trellis(6, [77 45]);

% Synchronization
Nper = M * Nsc/2;
metric = zeros(length(yFiltered) - 2*Nper + 1, 1);
for d = 1:length(metric)
    a = yFiltered(d:d+Nper-1);
    b = yFiltered(d+Nper:d+2*Nper-1);

    metric(d) = abs(sum(a .* conj(b))) / ...
        sqrt(sum(abs(a).^2) * sum(abs(b).^2));
end

[peakValue, peakIndex] = max(metric);

threshold = 0.8 * peakValue;

left = peakIndex;
while left > 1 && metric(left) > threshold
    left = left - 1;
end

right = peakIndex;
while right < length(metric) && metric(right) > threshold
    right = right + 1;
end

plateauStart = left + 1;

guard = round(M * Ncp / 2);
startSample = plateauStart + guard;

y = yFiltered(startSample:M:end);

% Decode
result = decode_ofdm(y, Nsc, Lsym, pilotTx, pilotIndex, trellis);

fprintf('\nDecoded message:\n');
disp(result.message);

% Plot synchronization metric
figure;
plot(metric);
grid on;
xlabel('Sample index before downsampling', 'FontSize', 18);
ylabel('Normalized periodic correlation', 'FontSize', 18);
title('Periodic Pilot Synchronization Metric', 'FontSize', 18);
set(gca, 'FontSize', 18)

% Plot constellation after equalization
figure;
plot(real(result.equalizedSymbols), imag(result.equalizedSymbols), '.');
grid on;
xlabel('In-phase', 'FontSize', 18);
ylabel('Quadrature', 'FontSize', 18);
title('Equalized QPSK Symbols', 'FontSize', 18);
set(gca, 'FontSize', 18)


function out = decode_ofdm(y, Nsc, Lsym, pilotTx, pilotIndex, trellis)

    out = struct();
    out.valid = false;
    out.message = '';
    out.equalizedSymbols = [];

    numBlocks = floor((length(y) - Nsc + 1) / Lsym);

    Yfreq = zeros(Nsc, numBlocks);

    for k = 1:numBlocks
        fftStart = 1 + (k-1)*Lsym;

        if fftStart + Nsc - 1 > length(y)
            break;
        end

        Yfreq(:, k) = fft(y(fftStart:fftStart+Nsc-1), Nsc);
    end

    % Channel estimation using the first OFDM pilot block
    Hpilot = Yfreq(pilotIndex, 1) ./ pilotTx(pilotIndex);

    % Interpolate missing subchannels
    Hhat = interp1(pilotIndex(:), Hpilot(:), (1:Nsc).', 'linear', 'extrap');

    % Equalization
    X = Yfreq ./ Hhat;

    % Decode message length from second OFDM block
    lengthSymbols = X(1:15, 2);
    lengthCodedBits = qpsk_demap(lengthSymbols);

    decodedLength = vitdec(lengthCodedBits, trellis, 15, 'term', 'hard');

    lengthBits = decodedLength(1:10);

    messageLength = lengthBits(:).' * (2.^(0:9)).';

    fprintf('Estimated message length = %d characters\n', messageLength);

    % Decode payload from 3rd OFDM block
    numMessageBits = 7 * messageLength;

    % Constraint length is 6, so encoder + 5 zeros
    numInputBits = numMessageBits + 5;

    payload = X(:, 3:end);
    payloadSymbols = payload(:);

    payloadSymbols = payloadSymbols(1:numInputBits);

    payloadCodedBits = qpsk_demap(payloadSymbols);

    decodedPayload = vitdec(payloadCodedBits, trellis, min(30, numInputBits), 'term', 'hard');

    decodedBits = decodedPayload(1:numMessageBits);

    out.valid = true;
    out.message = bits_to_ascii(decodedBits, messageLength);
    out.equalizedSymbols = payloadSymbols;
end


function bits = qpsk_demap(symbols)

    symbols = symbols(:);

    bits = zeros(1, 2*length(symbols));
    
    bits(1:2:end) = real(symbols) < 0;
    bits(2:2:end) = imag(symbols) < 0;
end


function message = bits_to_ascii(bits, messageLength)

    weights = 2.^(6:-1:0);

    message = char(zeros(1, messageLength));

    for k = 1:messageLength
        message(k) = char(bits(7*(k-1)+1:7*k) * weights.');
    end
end