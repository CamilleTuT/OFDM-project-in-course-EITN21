clear; clc; close all;

Fs = 44100;
fc = 10000;
Nsc = 128;
Ncp = 20;
Tsym = 58e-3;
M = round(Fs * Tsym / Nsc);
Lsym = Nsc + Ncp;

signalFile = sprintf('C:/Users/Hogwa/Desktop/SWE/4 PWC/Project_Part1/Signals_task3/signal1.mat');

data = load(signalFile);
r = data.R(:);
t = (0:length(r)-1).' / Fs;

baseband = 2 * r .* exp(-1i * 2*pi*fc*t);
[B, A] = butter(8, 0.05);
yFiltered = filter(B, A, baseband);

randn('state', 100);
P = sign(randn(1, Nsc/2));
pilotTx = zeros(Nsc, 1);
pilotTx(1:2:end) = 2 * P(:);

pilotIndex = 1:2:Nsc;
zeroIndex = 2:2:Nsc;
trellis = poly2trellis(6, [77 45]);

[candidates, yByPhase] = find_sync_candidates(yFiltered, M, Nsc, Ncp, pilotIndex, zeroIndex);
best = [];

for k = 1:length(candidates)
    cand = candidates(k);
    y = yByPhase{cand.phase};
    startRange = max(1, cand.start-Ncp-4) : min(length(y)-Nsc+1, cand.start+4);

    for fftStart = startRange
        trial = decode_from_start(y, fftStart, Nsc, Lsym, pilotTx, pilotIndex, trellis);
        if trial.valid && (isempty(best) || trial.metric < best.metric)
            best = trial;
        end
    end
end

fprintf('Decoded message:\n');
disp(best.message);

function [candidates, yByPhase] = find_sync_candidates(x, M, Nsc, Ncp, pilotIndex, zeroIndex)
    candidates = struct('score', {}, 'phase', {}, 'start', {});
    yByPhase = cell(M, 1);

    for phase = 1:M
        y = x(phase:M:end);
        yByPhase{phase} = y;
        maxStart = length(y) - Nsc + 1;
        if maxStart < 1
            continue;
        end

        energy = conv(abs(y).^2, ones(Nsc, 1), 'valid');
        energyThreshold = 0.05 * max(energy);
        scores = -inf(maxStart, 1);

        for startIndex = 1:maxStart
            if energy(startIndex) < energyThreshold
                continue;
            end

            Y = fft(y(startIndex:startIndex+Nsc-1), Nsc);
            pilotScore = mean(abs(Y(pilotIndex)).^2) / (mean(abs(Y(zeroIndex)).^2) + 1e-12);

            cpScore = 0;
            cpStart = startIndex - Ncp;
            if cpStart >= 1
                cp = y(cpStart:cpStart+Ncp-1);
                tail = y(startIndex+Nsc-Ncp:startIndex+Nsc-1);
                c = sum(conj(cp) .* tail);
                e = sum(abs(cp).^2) + sum(abs(tail).^2) + 1e-12;
                cpScore = abs(c)^2 / e^2;
            end

            scores(startIndex) = pilotScore + 5*cpScore - 1e-4*startIndex;
        end

        validScores = scores(isfinite(scores));
        if isempty(validScores)
            continue;
        end

        idx = find(scores >= 0.90 * max(validScores));
        for n = 1:length(idx)
            candidates(end+1).score = scores(idx(n));
            candidates(end).phase = phase;
            candidates(end).start = idx(n);
        end
    end

    [~, order] = sort([candidates.score], 'descend');
    candidates = candidates(order(1:min(80, length(order))));
end

function out = decode_from_start(y, fftStart, Nsc, Lsym, pilotTx, pilotIndex, trellis)
    out = struct('valid', false, 'metric', inf, 'message', '');
    numBlocks = floor((length(y) - fftStart + 1) / Lsym);
    if numBlocks < 3
        return;
    end

    Yfreq = zeros(Nsc, numBlocks);
    for k = 1:numBlocks
        idx = fftStart + (k-1)*Lsym;
        Yfreq(:, k) = fft(y(idx:idx+Nsc-1), Nsc);
    end

    Hpilot = Yfreq(pilotIndex, 1) ./ pilotTx(pilotIndex);
    Hhat = interp1(pilotIndex(:), Hpilot(:), (1:Nsc).', 'linear', 'extrap');
    if any(~isfinite(Hhat)) || any(abs(Hhat) < 1e-12)
        return;
    end

    X = Yfreq ./ Hhat;

    lengthSymbols = X(1:15, 2);
    lengthCodedBits = qpsk_demap(lengthSymbols);
    decodedLength = vitdec(lengthCodedBits, trellis, 15, 'term', 'hard');
    lengthBits = decodedLength(1:10);
    messageLength = lengthBits(:).' * (2.^(0:9)).';

    if messageLength < 1 || messageLength > 200
        return;
    end

    numMessageBits = 7 * messageLength;
    numInputBits = numMessageBits + 5;
    payload = X(:, 3:end);
    payloadSymbols = payload(:);
    if numInputBits > length(payloadSymbols)
        return;
    end

    payloadSymbols = payloadSymbols(1:numInputBits);
    payloadCodedBits = qpsk_demap(payloadSymbols);
    decodedPayload = vitdec(payloadCodedBits, trellis, min(30, numInputBits), 'term', 'hard');
    decodedBits = decodedPayload(1:numMessageBits);

    lengthErrors = sum(convenc([lengthBits(:).' zeros(1, 5)], trellis) ~= lengthCodedBits);
    payloadErrors = sum(convenc(decodedPayload(:).', trellis) ~= payloadCodedBits);

    out.valid = true;
    out.metric = 1000*lengthErrors + payloadErrors + 20*qpsk_distance([lengthSymbols; payloadSymbols]);
    out.message = bits_to_ascii(decodedBits, messageLength);
end

function bits = qpsk_demap(symbols)
    symbols = symbols(:);
    bits = zeros(1, 2*length(symbols));
    bits(1:2:end) = real(symbols) < 0;
    bits(2:2:end) = imag(symbols) < 0;
end

function d = qpsk_distance(symbols)
    nearest = (1 - 2*(real(symbols) < 0)) + 1i*(1 - 2*(imag(symbols) < 0));
    d = mean(abs(symbols(:) - nearest(:)).^2);
end

function message = bits_to_ascii(bits, messageLength)
    weights = 2.^(6:-1:0);
    message = char(zeros(1, messageLength));
    for k = 1:messageLength
        message(k) = char(bits(7*(k-1)+1:7*k) * weights.');
    end
end
