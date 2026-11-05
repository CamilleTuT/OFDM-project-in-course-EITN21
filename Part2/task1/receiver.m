clear; clc; close all;

% ===================== User settings =====================
recordedFile = 'recorded_task1_audio.wav';
outputFile = 'recovered_task1_text.txt';
txMatFile = 'task1_tx.mat';
maxSyncCandidates = 30;
figureFolder = 'task1_figures_rx';

% ===================== OFDM/audio parameters =====================
Fs = 44100;
fc = 10000;
Nfft = 128;
Ncp = 20;
M = 20;

activeBins = [(-50:-1), (1:50)].';
pilotBins = [-33; -11; 11; 33];
activeIndex = activeBins + Nfft/2 + 1;
pilotIndex = pilotBins + Nfft/2 + 1;
dataIndex = setdiff(activeIndex, pilotIndex, 'stable');
pilotSymbols = [1+1i; 1-1i; -1+1i; -1-1i];

rxLowpassCutoff = 0.055;
preambleSeed = 100;

constraintLength = 6;
codeGenerators = [77 45];
viterbiTraceback = 30;
headerInfoBits = 32;

cfg = build_task1_config(Fs, fc, Nfft, Ncp, M, ...
    activeIndex, pilotIndex, dataIndex, pilotSymbols, ...
    rxLowpassCutoff, preambleSeed, constraintLength, ...
    codeGenerators, viterbiTraceback, headerInfoBits);

cfg.recordStartTime = datetime(2026, 6, 12, 18, 03, 30);
cfg.recordSeconds = 35;

trellis = poly2trellis(cfg.constraintLength, cfg.codeGenerators);

wait_until_record_start(cfg.recordStartTime);

recorder = audiorecorder(cfg.Fs, 16, 1);
fprintf('Recording for %.1f seconds...\n', cfg.recordSeconds);
recordblocking(recorder, cfg.recordSeconds);
rxAudio = getaudiodata(recorder, 'double');
audiowrite(recordedFile, rxAudio, cfg.Fs);

make_receiver_input_plots(rxAudio, cfg, figureFolder);

[candidates, yByPhase] = find_sync_candidates_audio( ...
    rxAudio, cfg, maxSyncCandidates);
if isempty(candidates)
    error('No repeated OFDM preamble was found.');
end

best = [];
for k = 1:numel(candidates)
    cand = candidates(k);
    trial = decode_candidate(yByPhase{cand.phase}, cand, cfg, trellis);
    if trial.valid && (isempty(best) || trial.metric < best.metric)
        best = trial;
    end
end

if isempty(best)
    error('Preamble candidates were found, but no valid frame decoded.');
end

receivedBytes = bits_to_bytes(best.payloadBits);
receivedBytes = receivedBytes(1:best.payloadByteCount);
write_file_bytes(outputFile, receivedBytes);

fprintf('\nReceiver finished.\n');
fprintf('Recovered file: %s\n', outputFile);
fprintf('Repeated-preamble metric: %.4f\n', best.repeatScore);
fprintf('Known-preamble correlation: %.4f\n', best.templateScore);
fprintf('Estimated CFO: %.4f Hz\n', best.cfoHz);
fprintf('Payload bytes: %d\n', best.payloadByteCount);

if exist(txMatFile, 'file')
    tx = load(txMatFile);
    print_receiver_comparison(best, receivedBytes, tx, cfg, trellis);
end

plot_constellation(best.constellation, figureFolder, ...
    'rx_constellation.png');


function wait_until_record_start(recordStartTime)
    currentTime = datetime('now');

    if recordStartTime <= currentTime
        fprintf('Configured start time %s has passed. Recording starts now.\n', ...
            char(recordStartTime));
        return;
    end

    fprintf('Waiting until %s to start recording...\n', ...
        char(recordStartTime));

    while true
        secondsLeft = seconds(recordStartTime - datetime('now'));
        if secondsLeft <= 0
            break;
        end
        pause(min(secondsLeft, 0.2));
    end
end


function cfg = build_task1_config(Fs, fc, Nfft, Ncp, M, ...
        activeIndex, pilotIndex, dataIndex, pilotSymbols, ...
        rxLowpassCutoff, preambleSeed, constraintLength, ...
        codeGenerators, viterbiTraceback, headerInfoBits)

    cfg.Fs = Fs;
    cfg.fc = fc;
    cfg.Nsc = Nfft;
    cfg.Ncp = Ncp;
    cfg.M = M;
    cfg.FsBaseband = cfg.Fs / cfg.M;
    cfg.Lsym = cfg.Nsc + cfg.Ncp;
    cfg.preambleLengthSamples = cfg.Ncp + 2*cfg.Nsc;

    cfg.constraintLength = constraintLength;
    cfg.tailBits = cfg.constraintLength - 1;
    cfg.codeGenerators = codeGenerators;
    cfg.viterbiTraceback = viterbiTraceback;
    cfg.headerInfoBits = headerInfoBits;

    cfg.rxLpfNormalizedCutoff = rxLowpassCutoff;
    cfg.activeShiftIndices = activeIndex(:);
    cfg.pilotShiftIndices = pilotIndex(:);
    cfg.dataShiftIndices = dataIndex(:);
    cfg.numDataCarriers = numel(cfg.dataShiftIndices);
    cfg.dataPilotValues = pilotSymbols(:);
    cfg.preambleSeed = preambleSeed;
end

function [candidates, yByPhase] = find_sync_candidates_audio( ...
        rxAudio, cfg, maxCandidates)

    rxAudio = rxAudio(:);
    t = (0:numel(rxAudio)-1).' / cfg.Fs;
    complexBaseband = 2 * rxAudio .* exp(-1i * 2*pi*cfg.fc*t);

    [B, A] = butter(8, cfg.rxLpfNormalizedCutoff);
    complexBaseband = filter(B, A, complexBaseband);

    yByPhase = cell(cfg.M, 1);
    for phase = 1:cfg.M
        yByPhase{phase} = complexBaseband(phase:cfg.M:end);
    end

    candidates = find_repeated_preamble_candidates( ...
        yByPhase, cfg, maxCandidates);
end

function candidates = find_repeated_preamble_candidates( ...
        yByPhase, cfg, maxCandidates)

    candidates = struct('score', {}, 'repeatScore', {}, ...
        'templateScore', {}, 'phase', {}, 'cpStart', {}, ...
        'cfoHz', {});

    N = cfg.Nsc;
    averagingWindow = ones(N, 1);
    template = make_preamble_time(cfg);
    templateEnergy = sum(abs(template).^2);

    for phase = 1:numel(yByPhase)
        y = yByPhase{phase};
        if numel(y) < cfg.preambleLengthSamples
            continue;
        end

        first = y(1:end-N);
        second = y(N+1:end);
        product = conj(first) .* second;

        correlation = conv(product, averagingWindow, 'valid');
        firstEnergy = conv(abs(first).^2, averagingWindow, 'valid');
        secondEnergy = conv(abs(second).^2, averagingWindow, 'valid');
        repeatMetric = abs(correlation).^2 ./ ...
            (firstEnergy .* secondEnergy + eps);

        [sortedMetric, order] = sort(repeatMetric, 'descend');
        keptStarts = [];
        keptForPhase = 0;

        for n = 1:numel(order)
            if keptForPhase >= maxCandidates
                break;
            end
            if sortedMetric(n) < 0.10 * sortedMetric(1)
                break;
            end

            repeatStart = order(n);
            if any(abs(repeatStart - keptStarts) < cfg.Ncp)
                continue;
            end

            initialCfo = angle(correlation(repeatStart)) * ...
                cfg.FsBaseband / (2*pi*N);
            [cpStart, templateScore] = refine_preamble_start( ...
                y, repeatStart, initialCfo, template, ...
                templateEnergy, cfg);

            if isempty(cpStart)
                continue;
            end

            p1Start = cpStart + cfg.Ncp;
            p1 = y(p1Start:p1Start+N-1);
            p2 = y(p1Start+N:p1Start+2*N-1);
            repeatedCorrelation = sum(conj(p1) .* p2);
            cfoHz = angle(repeatedCorrelation) * ...
                cfg.FsBaseband / (2*pi*N);

            repeatScore = abs(repeatedCorrelation)^2 / ...
                ((sum(abs(p1).^2) * sum(abs(p2).^2)) + eps);

            candidates(end+1).score = repeatScore + templateScore; %#ok<AGROW>
            candidates(end).repeatScore = repeatScore;
            candidates(end).templateScore = templateScore;
            candidates(end).phase = phase;
            candidates(end).cpStart = cpStart;
            candidates(end).cfoHz = cfoHz;

            keptStarts(end+1) = repeatStart; %#ok<AGROW>
            keptForPhase = keptForPhase + 1;
        end
    end

    if ~isempty(candidates)
        [~, order] = sort([candidates.score], 'descend');
        candidates = candidates(order(1:min(maxCandidates, numel(order))));
    end
end

function [bestStart, bestScore] = refine_preamble_start( ...
        y, repeatStart, cfoHz, template, templateEnergy, cfg)

    searchFirst = max(1, repeatStart-cfg.Ncp-6);
    searchLast = min(repeatStart+6, ...
        numel(y)-cfg.preambleLengthSamples+1);

    bestStart = [];
    bestScore = -inf;

    for cpStart = searchFirst:searchLast
        idx = (cpStart:cpStart+cfg.preambleLengthSamples-1).';
        corrected = y(idx) .* exp(-1i * 2*pi*cfoHz * ...
            (idx-1) / cfg.FsBaseband);

        score = abs(sum(conj(template) .* corrected))^2 / ...
            ((sum(abs(corrected).^2) + eps) * templateEnergy);

        if score > bestScore
            bestScore = score;
            bestStart = cpStart;
        end
    end
end

function out = decode_candidate(y, candidate, cfg, trellis)
    out = struct('valid', false, 'metric', inf, ...
        'payloadBits', [], 'payloadByteCount', 0, ...
        'syncScore', candidate.score, 'cfoHz', candidate.cfoHz, ...
        'repeatScore', candidate.repeatScore, ...
        'templateScore', candidate.templateScore, ...
        'constellation', [], 'payloadCodedBits', []);

    sampleNumber = (0:numel(y)-1).';
    y = y .* exp(-1i * 2*pi*candidate.cfoHz * ...
        sampleNumber / cfg.FsBaseband);

    p1Start = candidate.cpStart + cfg.Ncp;
    p2Start = p1Start + cfg.Nsc;
    dataFftStart = p2Start + cfg.Nsc + cfg.Ncp;

    if dataFftStart + cfg.Nsc - 1 > numel(y)
        return;
    end

    preambleShift = make_preamble_shift(cfg);
    Y1 = read_useful_symbol(y, p1Start, cfg);
    Y2 = read_useful_symbol(y, p2Start, cfg);
    H1 = Y1(cfg.activeShiftIndices) ./ ...
        preambleShift(cfg.activeShiftIndices);
    H2 = Y2(cfg.activeShiftIndices) ./ ...
        preambleShift(cfg.activeShiftIndices);
    Hactive = (H1 + H2) / 2;

    blocksAvailable = 1 + floor( ...
        (numel(y) - (dataFftStart + cfg.Nsc - 1)) / cfg.Lsym);
    if blocksAvailable < 2
        return;
    end

    headerQpskCount = cfg.headerInfoBits + cfg.tailBits;
    [headerSymbols, headerDist] = extract_data_symbols( ...
        y, dataFftStart, 1, 1, headerQpskCount, cfg, Hactive);

    headerCodedBits = qpsk_demap(headerSymbols);
    headerDecoded = vitdec(headerCodedBits, trellis, ...
        min(cfg.viterbiTraceback, numel(headerSymbols)), ...
        'term', 'hard');
    headerBits = headerDecoded(1:cfg.headerInfoBits);
    payloadByteCount = bits_to_uint(headerBits);

    if payloadByteCount < 1
        return;
    end

    payloadBitCount = 8 * payloadByteCount;
    payloadQpskCount = payloadBitCount + cfg.tailBits;
    payloadBlocks = ceil(payloadQpskCount / cfg.numDataCarriers);

    if blocksAvailable < 1 + payloadBlocks
        return;
    end

    [payloadSymbols, payloadDist] = extract_data_symbols( ...
        y, dataFftStart, 2, payloadBlocks, payloadQpskCount, ...
        cfg, Hactive);

    payloadCodedBits = qpsk_demap(payloadSymbols);
    payloadDecoded = vitdec(payloadCodedBits, trellis, ...
        min(cfg.viterbiTraceback, numel(payloadSymbols)), ...
        'term', 'hard');
    payloadBits = payloadDecoded(1:payloadBitCount);

    headerReencoded = convenc(headerDecoded(:).', trellis);
    payloadReencoded = convenc(payloadDecoded(:).', trellis);
    headerErrors = sum(headerReencoded ~= headerCodedBits);
    payloadErrors = sum(payloadReencoded ~= payloadCodedBits);

    out.valid = true;
    out.metric = 1000*headerErrors + payloadErrors + ...
        20*(headerDist + payloadDist) - 10*candidate.score;
    out.payloadBits = payloadBits;
    out.payloadByteCount = payloadByteCount;
    out.payloadCodedBits = payloadCodedBits;
    out.constellation = [headerSymbols(:); ...
        payloadSymbols(1:min(end, 2000))];
end

function Yshift = read_useful_symbol(y, startIndex, cfg)
    noCp = y(startIndex:startIndex+cfg.Nsc-1);
    Yshift = fftshift(fft(noCp, cfg.Nsc)) / sqrt(cfg.Nsc);
end

function [symbols, distance] = extract_data_symbols( ...
        y, dataFftStart, firstBlock, numBlocks, needed, cfg, Hactive)

    symbols = zeros(numBlocks * cfg.numDataCarriers, 1);
    allBins = cfg.activeShiftIndices - (cfg.Nsc/2 + 1);
    pilotBins = cfg.pilotShiftIndices - (cfg.Nsc/2 + 1);
    pos = 1;

    for b = 1:numBlocks
        blockNumber = firstBlock + b - 1;
        symbolStart = dataFftStart + (blockNumber-1)*cfg.Lsym;
        Yshift = read_useful_symbol(y, symbolStart, cfg);

        Xeq = zeros(cfg.Nsc, 1);
        Xeq(cfg.activeShiftIndices) = ...
            Yshift(cfg.activeShiftIndices) ./ Hactive;

        pilotGain = Xeq(cfg.pilotShiftIndices) ./ ...
            cfg.dataPilotValues;
        pilotPhase = unwrap(angle(pilotGain));
        phaseFit = polyfit(double(pilotBins), pilotPhase, 1);
        phaseCorrection = exp(-1i * ...
            polyval(phaseFit, double(allBins)));
        Xeq(cfg.activeShiftIndices) = ...
            Xeq(cfg.activeShiftIndices) .* phaseCorrection;

        symbols(pos:pos+cfg.numDataCarriers-1) = ...
            Xeq(cfg.dataShiftIndices);
        pos = pos + cfg.numDataCarriers;
    end

    symbols = symbols(1:needed);
    distance = qpsk_distance(symbols);
end

function bits = qpsk_demap(symbols)
    symbols = symbols(:);
    bits = zeros(1, 2*numel(symbols));
    bits(1:2:end) = real(symbols) < 0;
    bits(2:2:end) = imag(symbols) < 0;
end

function distance = qpsk_distance(symbols)
    nearest = (1 - 2*(real(symbols) < 0)) + ...
        1i*(1 - 2*(imag(symbols) < 0));
    distance = mean(abs(symbols(:) - nearest(:)).^2);
end

function bytes = bits_to_bytes(bits)
    bits = bits(:).';
    weights = 2.^(7:-1:0);
    bytes = zeros(1, numel(bits)/8, 'uint8');
    for k = 1:numel(bytes)
        bytes(k) = uint8(bits(8*(k-1)+1:8*k) * weights.');
    end
end

function value = bits_to_uint(bits)
    value = uint64(0);
    for k = 1:numel(bits)
        value = bitshift(value, 1) + uint64(bits(k) ~= 0);
    end
    value = double(value);
end

function write_file_bytes(fileName, bytes)
    fid = fopen(fileName, 'wb');
    fwrite(fid, uint8(bytes), 'uint8');
    fclose(fid);
end

function Xshift = make_preamble_shift(cfg)
    rng(cfg.preambleSeed);
    Xshift = zeros(cfg.Nsc, 1);
    knownBpsk = 1 - 2*(rand( ...
        numel(cfg.activeShiftIndices), 1) > 0.5);
    Xshift(cfg.activeShiftIndices) = sqrt(2) * knownBpsk;
end

function preamble = make_preamble_time(cfg)
    preambleShift = make_preamble_shift(cfg);
    pilotOfdmSymbol = sqrt(cfg.Nsc) * ...
        ifft(ifftshift(preambleShift), cfg.Nsc);
    preambleCp = pilotOfdmSymbol(end-cfg.Ncp+1:end);
    preamble = [preambleCp; pilotOfdmSymbol; pilotOfdmSymbol];
end

function make_receiver_input_plots(rxAudio, cfg, figureFolder)
    if ~exist(figureFolder, 'dir')
        mkdir(figureFolder);
    end
    plot_time_signal(rxAudio, cfg.Fs, ...
        'Received audio waveform', 'Time / s', 'Amplitude', ...
        figureFolder, 'rx_audio_waveform.png');
    plot_spectrum(rxAudio, cfg.Fs, ...
        'Received audio spectrum', figureFolder, ...
        'rx_audio_spectrum.png');
end

function print_receiver_comparison(best, receivedBytes, tx, cfg, trellis)
    if isfield(tx, 'sourceBits') && ~isempty(best.payloadCodedBits)
        txCodedBits = convenc([tx.sourceBits(:).' ...
            zeros(1, cfg.tailBits)], trellis);
        n = min(numel(txCodedBits), numel(best.payloadCodedBits));
        bitErrors = sum(txCodedBits(1:n) ~= ...
            best.payloadCodedBits(1:n));
        fprintf('Raw coded BER before Viterbi: %.6f (%d / %d)\n', ...
            bitErrors/n, bitErrors, n);
    end

    if isfield(tx, 'payloadBytes')
        n = min(numel(tx.payloadBytes), numel(receivedBytes));
        byteErrors = sum(uint8(tx.payloadBytes(1:n)) ~= ...
            uint8(receivedBytes(1:n)));
        fprintf('Byte errors after decoding: %d / %d\n', ...
            byteErrors, n);
    end
end

function plot_time_signal(x, Fs, plotTitle, ...
        xLabelText, yLabelText, folder, fileName)

    t = (0:numel(x)-1).' / Fs;
    f = figure('Visible', 'off');
    plot(t, real(x));
    grid on;
    xlabel(xLabelText);
    ylabel(yLabelText);
    title(plotTitle);
    saveas(f, fullfile(folder, fileName));
    close(f);
end

function plot_spectrum(x, Fs, plotTitle, folder, fileName)
    Nplot = 2^nextpow2(min(numel(x), 262144));
    X = fftshift(fft(x, Nplot));
    fAxis = (-Nplot/2:Nplot/2-1).' * Fs / Nplot;
    magDb = 20*log10(abs(X)/max(abs(X)+eps) + eps);

    f = figure('Visible', 'off');
    plot(fAxis, magDb);
    grid on;
    xlabel('Frequency / Hz');
    ylabel('Magnitude / dB');
    title(plotTitle);
    ylim([-100 5]);
    saveas(f, fullfile(folder, fileName));
    close(f);
end

function plot_constellation(symbols, folder, fileName)
    if ~exist(folder, 'dir')
        mkdir(folder);
    end
    f = figure('Visible', 'off');
    plot(real(symbols), imag(symbols), '.');
    axis equal;
    grid on;
    xlabel('In-phase');
    ylabel('Quadrature');
    title('Equalized QPSK constellation');
    saveas(f, fullfile(folder, fileName));
    close(f);
end
