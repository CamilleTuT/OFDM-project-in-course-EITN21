clear; clc; close all;

outputFile = 'recovered_task1_text.txt';
txMatFile = 'task1_tx.mat';
figureFolder = 'task1_figures_rx';

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

cfg.recordStartTime = datetime(2026, 6, 12, 23, 24, 40);
cfg.recordSeconds = 35;

trellis = poly2trellis(cfg.constraintLength, cfg.codeGenerators);

wait_until_record_start(cfg.recordStartTime);

recorder = audiorecorder(cfg.Fs, 16, 1);
fprintf('Recording for %.1f seconds...\n', cfg.recordSeconds);
recordblocking(recorder, cfg.recordSeconds);
rxAudio = getaudiodata(recorder, 'double');

make_receiver_input_plots(rxAudio, cfg, figureFolder);

[y, yBeforeCfo, sync] = synchronize_full_rate(rxAudio, cfg);
best = decode_synchronized_frame(y, sync, cfg, trellis);

if ~best.valid
    error('The synchronized OFDM frame could not be decoded.');
end

receivedBytes = bits_to_bytes(best.payloadBits);
write_file_bytes(outputFile, receivedBytes);

fprintf('\nReceiver finished.\n');
fprintf('Correlation peak value: %.4f\n', best.peakValue);
fprintf('Estimated CFO: %.4f Hz\n', best.cfoHz);
fprintf('Correlation peak: audio sample %d\n', sync.peakIndex);
fprintf('Payload bytes: %d\n', best.payloadByteCount);

tx = load(txMatFile);
print_receiver_comparison(best, receivedBytes, tx, cfg, trellis);

plot_constellation(best.constellation, figureFolder, 'rx_constellation.png');
plot_pilot_phase_evolution(yBeforeCfo, best.pilotPhaseAfterCompensation, ...
    best.payloadByteCount, cfg, figureFolder);


function wait_until_record_start(recordStartTime)
    fprintf('Waiting until %s to start recording...\n', char(recordStartTime));
    secondsLeft = seconds(recordStartTime - datetime('now'));
    pause(secondsLeft);
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

function [y, yBeforeCfo, sync] = synchronize_full_rate(rxAudio, cfg)
    rxAudio = rxAudio(:);
    audioSample = (0:numel(rxAudio)-1).';
    % down conversion
    complexBaseband = 2 * rxAudio .* exp(-1i * 2*pi*cfg.fc * audioSample / cfg.Fs);
    [B, A] = butter(8, cfg.rxLpfNormalizedCutoff);
    complexBaseband = filter(B, A, complexBaseband);

    repeatDistance = cfg.M * cfg.Nsc;
    numPositions = length(complexBaseband) - 2*repeatDistance + 1;
    correlation = zeros(numPositions, 1);
    firstEnergy = zeros(numPositions, 1);
    secondEnergy = zeros(numPositions, 1);
    metric = zeros(numPositions, 1);

    for d = 1:numPositions
        first = complexBaseband(d:d+repeatDistance-1);
        second = complexBaseband(d+repeatDistance:d+2*repeatDistance-1);

        correlation(d) = sum(conj(first) .* second);
        firstEnergy(d) = sum(abs(first).^2);
        secondEnergy(d) = sum(abs(second).^2);

        metric(d) = abs(correlation(d))^2 / (firstEnergy(d) * secondEnergy(d) + eps);
    end
    [peakValue, peakIndex] = max(metric);
    startSample = peakIndex;

    cfoHz = angle(correlation(startSample)) * cfg.Fs / (2*pi*repeatDistance);

    yBeforeCfo = complexBaseband(startSample:cfg.M:end);
    complexBaseband = complexBaseband .* exp(-1i * 2*pi*cfoHz * audioSample / cfg.Fs);
    y = complexBaseband(startSample:cfg.M:end);

    sync.peakValue = peakValue;
    sync.cfoHz = cfoHz;
    sync.peakIndex = peakIndex;
end

function out = decode_synchronized_frame(y, sync, cfg, trellis)
    out = struct('valid', false, ...
        'payloadBits', [], 'payloadByteCount', 0, ...
        'cfoHz', sync.cfoHz, 'peakValue', sync.peakValue, ...
        'constellation', [], 'payloadCodedBits', [], ...
        'pilotPhaseAfterCompensation', []);

    p1Start = 1;
    p2Start = p1Start + cfg.Nsc;
    dataFftStart = p2Start + cfg.Nsc + cfg.Ncp;

    if dataFftStart + cfg.Nsc - 1 > numel(y)
        return;
    end

    preambleShift = make_preamble_shift(cfg);
    Y1 = read_useful_symbol(y, p1Start, cfg);
    Y2 = read_useful_symbol(y, p2Start, cfg);
    H1 = Y1(cfg.activeShiftIndices) ./ preambleShift(cfg.activeShiftIndices);
    H2 = Y2(cfg.activeShiftIndices) ./ preambleShift(cfg.activeShiftIndices);
    Hactive = (H1 + H2) / 2;

    blocksAvailable = 1 + floor((numel(y) - (dataFftStart + cfg.Nsc - 1)) / cfg.Lsym);
    if blocksAvailable < 2
        return;
    end

    headerQpskCount = cfg.headerInfoBits + cfg.tailBits;
    [headerSymbols, headerPilotPhase] = extract_data_symbols( ...
        y, dataFftStart, 1, 1, headerQpskCount, cfg, Hactive);

    headerCodedBits = qpsk_demap(headerSymbols);
    headerDecoded = vitdec(headerCodedBits, trellis, cfg.viterbiTraceback, 'term', 'hard');
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

    [payloadSymbols, payloadPilotPhase] = extract_data_symbols( ...
        y, dataFftStart, 2, payloadBlocks, payloadQpskCount, cfg, Hactive);

    payloadCodedBits = qpsk_demap(payloadSymbols);
    payloadDecoded = vitdec(payloadCodedBits, trellis, cfg.viterbiTraceback, 'term', 'hard');
    payloadBits = payloadDecoded(1:payloadBitCount);

    out.valid = true;
    out.payloadBits = payloadBits;
    out.payloadByteCount = payloadByteCount;
    out.payloadCodedBits = payloadCodedBits;
    out.constellation = [headerSymbols(:); payloadSymbols(1:min(end, 2000))];
    out.pilotPhaseAfterCompensation = [headerPilotPhase, payloadPilotPhase];
end

function Yshift = read_useful_symbol(y, startIndex, cfg)
    noCp = y(startIndex:startIndex+cfg.Nsc-1);
    Yshift = fftshift(fft(noCp, cfg.Nsc)) / sqrt(cfg.Nsc);
end

function [symbols, pilotPhaseAfterCompensation] = extract_data_symbols( ...
        y, dataFftStart, firstBlock, numBlocks, needed, cfg, Hactive)

    symbols = zeros(numBlocks * cfg.numDataCarriers, 1);
    pilotPhaseAfterCompensation = zeros(numel(cfg.pilotShiftIndices), numBlocks);
    allBins = cfg.activeShiftIndices - (cfg.Nsc/2 + 1);
    pilotBins = cfg.pilotShiftIndices - (cfg.Nsc/2 + 1);
    pos = 1;

    for b = 1:numBlocks
        blockNumber = firstBlock + b - 1;
        symbolStart = dataFftStart + (blockNumber-1)*cfg.Lsym;
        Yshift = read_useful_symbol(y, symbolStart, cfg);

        Xeq = zeros(cfg.Nsc, 1);
        Xeq(cfg.activeShiftIndices) = Yshift(cfg.activeShiftIndices) ./ Hactive;

        pilotGain = Xeq(cfg.pilotShiftIndices) ./ cfg.dataPilotValues;
        pilotPhase = unwrap(angle(pilotGain));
        phaseFit = polyfit(double(pilotBins), pilotPhase, 1);
        phaseCorrection = exp(-1i * polyval(phaseFit, double(allBins)));
        Xeq(cfg.activeShiftIndices) = Xeq(cfg.activeShiftIndices) .* phaseCorrection;

        compensatedPilots = Xeq(cfg.pilotShiftIndices) ./ cfg.dataPilotValues;
        pilotPhaseAfterCompensation(:, b) = angle(compensatedPilots);
        symbols(pos:pos+cfg.numDataCarriers-1) = Xeq(cfg.dataShiftIndices);
        pos = pos + cfg.numDataCarriers;
    end

    symbols = symbols(1:needed);
end

function bits = qpsk_demap(symbols)
    symbols = symbols(:);
    bits = zeros(1, 2*numel(symbols));
    bits(1:2:end) = real(symbols) < 0;
    bits(2:2:end) = imag(symbols) < 0;
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
    knownBpsk = 1 - 2*(rand(numel(cfg.activeShiftIndices), 1) > 0.5);
    Xshift(cfg.activeShiftIndices) = sqrt(2) * knownBpsk;
end

function make_receiver_input_plots(rxAudio, cfg, figureFolder)
    if ~exist(figureFolder, 'dir')
        mkdir(figureFolder);
    end

    t = (0:numel(rxAudio)-1).' / cfg.Fs;
    f = figure();
    plot(t, rxAudio);
    grid on;
    xlabel('Time / s', 'FontSize', 18);
    ylabel('Amplitude', 'FontSize', 18);
    title('Received audio waveform', 'FontSize', 18);
    set(gca, 'FontSize', 18)
    saveas(f, fullfile(figureFolder, 'rx_audio_waveform.png'));
    close(f);

    Nplot = 2^nextpow2(min(numel(rxAudio), 262144));
    spectrum = fftshift(fft(rxAudio, Nplot));
    fAxis = (-Nplot/2:Nplot/2-1).' * cfg.Fs / Nplot;
    magDb = 20*log10(abs(spectrum)/max(abs(spectrum)+eps) + eps);

    f = figure();
    plot(fAxis, magDb);
    grid on;
    xlabel('Frequency / Hz', 'FontSize', 18);
    ylabel('Magnitude / dB', 'FontSize', 18);
    title('Received audio spectrum', 'FontSize', 18);
    set(gca, 'FontSize', 18)
    ylim([-100 5]);
    saveas(f, fullfile(figureFolder, 'rx_audio_spectrum.png'));
    close(f);
end

function print_receiver_comparison(best, receivedBytes, tx, cfg, trellis)
    txCodedBits = convenc([tx.sourceBits(:).' zeros(1, cfg.tailBits)], trellis);
    n = min(numel(txCodedBits), numel(best.payloadCodedBits));
    bitErrors = sum(txCodedBits(1:n) ~= best.payloadCodedBits(1:n));
    fprintf('Raw coded BER before Viterbi: %.6f (%d / %d)\n', bitErrors/n, bitErrors, n);

    n = min(numel(tx.payloadBytes), numel(receivedBytes));
    byteErrors = sum(uint8(tx.payloadBytes(1:n)) ~= uint8(receivedBytes(1:n)));
    fprintf('Byte errors after decoding: %d / %d\n', byteErrors, n);
end

function plot_constellation(symbols, folder, fileName)
    if ~exist(folder, 'dir')
        mkdir(folder);
    end
    f = figure('Visible', 'off');
    plot(real(symbols), imag(symbols), '.');
    axis equal;
    grid on;
    xlabel('In-phase', 'FontSize', 18);
    ylabel('Quadrature', 'FontSize', 18);
    title('Equalized QPSK constellation', 'FontSize', 18);
    set(gca, 'FontSize', 18)
    saveas(f, fullfile(folder, fileName));
    close(f);
end

function plot_pilot_phase_evolution(yBeforeCfo, phaseAfterCompensation, ...
        payloadByteCount, cfg, folder)
    dataFftStart = 2*cfg.Nsc + cfg.Ncp + 1;
    payloadQpskCount = 8*payloadByteCount + cfg.tailBits;
    numBlocks = 1 + ceil(payloadQpskCount / cfg.numDataCarriers);
    numPilots = numel(cfg.pilotShiftIndices);

    phaseBeforeCfo = zeros(numPilots, numBlocks);

    for b = 1:numBlocks
        symbolStart = dataFftStart + (b-1)*cfg.Lsym;
        Ybefore = read_useful_symbol(yBeforeCfo, symbolStart, cfg);

        pilotBefore = Ybefore(cfg.pilotShiftIndices) ./ cfg.dataPilotValues;
        phaseBeforeCfo(:, b) = angle(pilotBefore);
    end

    phaseBeforeCfo = unwrap(phaseBeforeCfo, [], 2);
    phaseAfterCompensation = unwrap(phaseAfterCompensation, [], 2);
    ofdmSymbolNumber = 1:numBlocks;

    f = figure;
    plot(ofdmSymbolNumber, phaseBeforeCfo.', 'LineWidth', 1);
    grid on;
    xlabel('OFDM symbol number', 'FontSize', 18);
    ylabel('Pilot phase / rad', 'FontSize', 18);
    title('Pilot phase evolution before CFO compensation', 'FontSize', 18);
    legend('Pilot -33', 'Pilot -11', 'Pilot 11', 'Pilot 33', 'FontSize', 16);
    set(gca, 'FontSize', 18);
    saveas(f, fullfile(folder, 'pilot_phase_before_cfo.png'));
    close(f);

    f = figure;
    plot(ofdmSymbolNumber, phaseAfterCompensation.', 'LineWidth', 1);
    grid on;
    xlabel('OFDM symbol number', 'FontSize', 18);
    ylabel('Pilot phase / rad', 'FontSize', 18);
    title('Pilot phase evolution after compensation', 'FontSize', 18);
    legend('Pilot -33', 'Pilot -11', 'Pilot 11', 'Pilot 33', 'FontSize', 16);
    set(gca, 'FontSize', 18);
    saveas(f, fullfile(folder, 'pilot_phase_after_compensation.png'));
    close(f);
end
