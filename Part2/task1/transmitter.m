clear; clc; close all;

inputFile = 'C:\Users\Hogwa\Desktop\SWE\4 PWC\Project_Part2\test1.txt';
txMatFile = 'task1_tx.mat';
txWavFile = 'task1_tx_audio.wav';
figureFolder = 'task1_figures_tx';

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

guardSymbols = 2;
audioAmplitude = 0.80;
rxLowpassCutoff = 0.055;
preambleSeed = 100;

constraintLength = 6;
codeGenerators = [77 45];
viterbiTraceback = 30;
headerInfoBits = 32;

cfg = build_task1_config(Fs, fc, Nfft, Ncp, M, ...
    activeIndex, pilotIndex, dataIndex, pilotSymbols, ...
    guardSymbols, audioAmplitude, ...
    rxLowpassCutoff, preambleSeed, constraintLength, codeGenerators, ...
    viterbiTraceback, headerInfoBits);

payloadBytes = read_file_bytes(inputFile);
sourceBits = bytes_to_bits(payloadBytes);

trellis = poly2trellis(cfg.constraintLength, cfg.codeGenerators);

headerBits = make_header_bits(numel(payloadBytes), cfg);
headerSymbols = encode_bits_to_qpsk(headerBits, trellis, cfg.tailBits);
payloadSymbols = encode_bits_to_qpsk(sourceBits, trellis, cfg.tailBits);

txFrame = build_ofdm_frame(headerSymbols, payloadSymbols, cfg);

guard = zeros(cfg.guardSymbols * cfg.Lsym, 1);
txBaseband = [guard; txFrame; guard];

txInterpolated = sinc_interpolate(txBaseband, cfg.M);
t = (0:numel(txInterpolated)-1).' / cfg.Fs;
txAudioFrame = real(txInterpolated .* exp(1i * 2*pi*cfg.fc*t));
txAudioFrame = cfg.audioAmplitude * txAudioFrame / max(abs(txAudioFrame));

audioGuard = zeros(round(2.0 * cfg.Fs), 1);
txAudio = [audioGuard; txAudioFrame; audioGuard];

audiowrite(txWavFile, txAudio, cfg.Fs);
save(txMatFile, 'cfg', 'inputFile', 'txWavFile', 'payloadBytes', ...
    'sourceBits', 'headerBits', 'headerSymbols', 'payloadSymbols', ...
    'txBaseband', 'txAudioFrame', 'txAudio');

make_transmitter_report_plots(txAudio, txAudioFrame, txBaseband, payloadBytes, cfg, figureFolder);

fprintf('Transmitter finished');


function cfg = build_task1_config(Fs, fc, Nfft, Ncp, M, ...
    activeIndex, pilotIndex, dataIndex, pilotSymbols, ...
    guardSymbols, audioAmplitude, ...
    rxLowpassCutoff, preambleSeed, constraintLength, codeGenerators, ...
    viterbiTraceback, headerInfoBits)

    cfg.Fs = Fs;
    cfg.fc = fc;
    cfg.Nsc = Nfft;
    cfg.Ncp = Ncp;
    cfg.M = M;
    cfg.FsBaseband = cfg.Fs / cfg.M;
    cfg.Lsym = cfg.Nsc + cfg.Ncp;

    cfg.preambleSpacingSamples = cfg.Nsc;
    cfg.preambleLengthSamples = cfg.Ncp + 2 * cfg.Nsc;
    cfg.guardSymbols = guardSymbols;

    cfg.constraintLength = constraintLength;
    cfg.tailBits = cfg.constraintLength - 1;
    cfg.codeGenerators = codeGenerators;
    cfg.viterbiTraceback = viterbiTraceback;

    cfg.headerInfoBits = headerInfoBits;

    cfg.audioAmplitude = audioAmplitude;
    cfg.rxLpfNormalizedCutoff = rxLowpassCutoff;

    cfg.activeShiftIndices = activeIndex(:);
    cfg.pilotShiftIndices = pilotIndex(:);
    cfg.dataShiftIndices = dataIndex(:);
    cfg.numDataCarriers = numel(cfg.dataShiftIndices);
    cfg.dataPilotValues = pilotSymbols(:);

    cfg.preambleSeed = preambleSeed;
    cfg.paddingSeed = preambleSeed + 1;
end

function bytes = read_file_bytes(fileName)
    fid = fopen(fileName, 'rb');
    bytes = fread(fid, Inf, '*uint8').';
    fclose(fid);
end

function bits = bytes_to_bits(bytes)
    bytes = uint8(bytes(:));
    bits = zeros(1, 8*numel(bytes));
    for b = 1:8
        bits(b:8:end) = bitget(bytes, 9-b);
    end
end

function bits = make_header_bits(payloadByteCount, cfg)
    bits = uint_to_bits(payloadByteCount, cfg.headerInfoBits);
end

function bits = uint_to_bits(value, width)
    value = uint64(value);
    bits = zeros(1, width);
    for k = 1:width
        bits(k) = bitget(value, width-k+1);
    end
end

function symbols = encode_bits_to_qpsk(infoBits, trellis, tailBits)
    codedBits = convenc([infoBits(:).' zeros(1, tailBits)], trellis);
    symbols = qpsk_mod(codedBits);
end

function symbols = qpsk_mod(bits)
    bits = bits(:).';
    if mod(numel(bits), 2) ~= 0
        bits(end+1) = 0;
    end
    pairs = reshape(bits, 2, []).';
    symbols = (1 - 2*pairs(:, 1)) + 1i*(1 - 2*pairs(:, 2));
end

function txFrame = build_ofdm_frame(headerSymbols, payloadSymbols, cfg)
    preambleShift = make_preamble_shift(cfg);
    pilotOfdmSymbol = make_ofdm_useful_symbol(preambleShift, cfg);
    preambleCp = pilotOfdmSymbol(end-cfg.Ncp+1:end);
    preamble = [preambleCp; pilotOfdmSymbol; pilotOfdmSymbol];

    txFrame = [preamble; ...
        make_data_ofdm_blocks(headerSymbols, cfg); ...
        make_data_ofdm_blocks(payloadSymbols, cfg)];
end

function blocks = make_data_ofdm_blocks(qpskSymbols, cfg)
    numBlocks = ceil(numel(qpskSymbols) / cfg.numDataCarriers);
    padCount = numBlocks * cfg.numDataCarriers - numel(qpskSymbols);
    paddingSymbols = make_qpsk_padding(padCount, cfg.paddingSeed);
    qpskSymbols = [qpskSymbols(:); paddingSymbols];

    blocks = zeros(numBlocks * cfg.Lsym, 1);
    pos = 1;
    for b = 1:numBlocks
        Xshift = zeros(cfg.Nsc, 1);
        Xshift(cfg.pilotShiftIndices) = cfg.dataPilotValues;
        Xshift(cfg.dataShiftIndices) = qpskSymbols(pos:pos+cfg.numDataCarriers-1);
        blocks((b-1)*cfg.Lsym+1:b*cfg.Lsym) = make_ofdm_time_symbol(Xshift, cfg);
        pos = pos + cfg.numDataCarriers;
    end
end

function symbols = make_qpsk_padding(count, seed)
    if count == 0
        symbols = complex(zeros(0, 1));
        return;
    end

    stream = RandStream('mt19937ar', 'Seed', seed);
    bits = randi(stream, [0 1], count, 2);
    symbols = (1 - 2*bits(:, 1)) + 1i*(1 - 2*bits(:, 2));
end

function Xshift = make_preamble_shift(cfg)
    rng(cfg.preambleSeed);
    Xshift = zeros(cfg.Nsc, 1);
    knownBpsk = 1 - 2*(rand(numel(cfg.activeShiftIndices), 1) > 0.5);
    Xshift(cfg.activeShiftIndices) = sqrt(2) * knownBpsk;
end

function xcp = make_ofdm_time_symbol(Xshift, cfg)
    x = make_ofdm_useful_symbol(Xshift, cfg);
    xcp = [x(end-cfg.Ncp+1:end); x];
end

function x = make_ofdm_useful_symbol(Xshift, cfg)
    x = sqrt(cfg.Nsc) * ifft(ifftshift(Xshift), cfg.Nsc);
end

function y = sinc_interpolate(x, M)
    x = x(:);
    up = zeros(numel(x) * M, 1);
    up(1:M:end) = x;

    span = 10;
    n = (-span*M:span*M).';
    h = local_sinc(n / M) .* hamming(numel(n));
    h = h / h(ceil(numel(h)/2));

    y = conv(up, h, 'same');
end

function y = local_sinc(x)
    y = ones(size(x));
    idx = abs(x) > 1e-12;
    y(idx) = sin(pi*x(idx)) ./ (pi*x(idx));
end


function make_transmitter_report_plots(txAudio, txAudioFrame, txBaseband, payloadBytes, cfg, figureFolder)
    if ~exist(figureFolder, 'dir')
        mkdir(figureFolder);
    end

    print_transmitter_numbers(txAudioFrame, payloadBytes, cfg);
    plot_time_signal(txAudio, cfg.Fs, 'Transmit audio waveform', 'Time / s', 'Amplitude', figureFolder, 'tx_audio_waveform.png');
    plot_spectrum(txAudio, cfg.Fs, 'Transmit audio spectrum', figureFolder, 'tx_audio_spectrum.png');
    plot_time_magnitude(txBaseband, cfg.FsBaseband, 'Complex baseband OFDM waveform', figureFolder, 'tx_baseband_waveform.png');
    plot_spectrum(txBaseband, cfg.FsBaseband, 'Complex baseband OFDM spectrum', figureFolder, 'tx_baseband_spectrum.png');
    plot_subcarrier_allocation(cfg, figureFolder, 'subcarrier_allocation.png');
end

function print_transmitter_numbers(txAudioFrame, payloadBytes, cfg)
    symbolRate = cfg.FsBaseband / cfg.Lsym;
    theoreticalCodedBitRate = symbolRate * cfg.numDataCarriers * 2;
    theoreticalInfoBitRate = theoreticalCodedBitRate / 2;
    transmissionDuration = numel(txAudioFrame) / cfg.Fs;
    practicalBitRate = 8 * numel(payloadBytes) / transmissionDuration;
    paprDb = 10*log10(max(abs(txAudioFrame).^2) / mean(abs(txAudioFrame).^2));

    fprintf('Baseband OFDM symbol rate: %.2f symbols/s\n', symbolRate);
    fprintf('Theoretical coded bit rate: %.1f bit/s\n', theoreticalCodedBitRate);
    fprintf('Theoretical information bit rate after 1/2 rate code: %.1f bit/s\n', theoreticalInfoBitRate);
    fprintf('Practical payload bit rate: %.1f bit/s\n', practicalBitRate);
    fprintf('Transmit PAPR: %.2f dB\n\n', paprDb);
end

function plot_time_signal(x, Fs, plotTitle, xLabelText, yLabelText, folder, fileName)
    t = (0:numel(x)-1).' / Fs;
    f = figure;
    plot(t, real(x));
    grid on;
    xlabel(xLabelText, 'FontSize', 18);
    ylabel(yLabelText, 'FontSize', 18);
    title(plotTitle, 'FontSize', 18);
    set(gca, 'FontSize', 18)
    saveas(f, fullfile(folder, fileName));
end

function plot_time_magnitude(x, Fs, plotTitle, folder, fileName)
    t = (0:numel(x)-1).' / Fs;
    f = figure;
    plot(t, abs(x));
    grid on;
    xlabel('Time / s', 'FontSize', 18);
    ylabel('Magnitude', 'FontSize', 18);
    title(plotTitle, 'FontSize', 18);
    set(gca, 'FontSize', 18)
    saveas(f, fullfile(folder, fileName));
end

function plot_spectrum(x, Fs, plotTitle, folder, fileName)
    Nplot = 2^nextpow2(min(numel(x), 262144));
    X = fftshift(fft(x, Nplot));
    fAxis = (-Nplot/2:Nplot/2-1).' * Fs / Nplot;
    magDb = 20*log10(abs(X)/max(abs(X)));

    f = figure;
    plot(fAxis, magDb);
    grid on;
    xlabel('Frequency / Hz', 'FontSize', 18);
    ylabel('Magnitude / dB', 'FontSize', 18);
    title(plotTitle, 'FontSize', 18);
    set(gca, 'FontSize', 18)
    ylim([-100 5]);
    saveas(f, fullfile(folder, fileName));
end

function plot_subcarrier_allocation(cfg, folder, fileName)
    bins = (-cfg.Nsc/2:cfg.Nsc/2-1).';
    values = zeros(cfg.Nsc, 1);
    values(cfg.activeShiftIndices) = 1;
    values(cfg.pilotShiftIndices) = 2;

    f = figure;
    stem(bins, values, 'filled');
    grid on;
    xlabel('Subcarrier index', 'FontSize', 18);
    ylabel('Allocation', 'FontSize', 18);
    title('OFDM subcarrier allocation', 'FontSize', 18);
    set(gca, 'FontSize', 18)
    ylim([-0.2 2.4]);
    saveas(f, fullfile(folder, fileName));
end
