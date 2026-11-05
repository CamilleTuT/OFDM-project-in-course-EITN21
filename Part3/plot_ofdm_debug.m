function plot_ofdm_debug(debug, cfg, mode)
%PLOT_OFDM_DEBUG Draw waveform, spectrum, sync, constellation, and pilot plots.

if ~cfg.enablePlots
    return;
end

mode = lower(char(mode));
switch mode
    case 'tx'
        if isfield(debug, 'txWaveform')
            figure('Name', 'TX time-domain waveform');
            plot(real(debug.txWaveform));
            grid on;
            xlabel('Sample');
            ylabel('Real amplitude');
            title('Transmit time-domain waveform');

            plot_spectrum(debug.txWaveform, cfg.sampleRate, 'Transmit spectrum');
        end

    case 'rx'
        if isfield(debug, 'rxSignal')
            figure('Name', 'RX time-domain waveform');
            plot(abs(debug.rxSignal));
            grid on;
            xlabel('Sample');
            ylabel('Magnitude');
            title('Received time-domain waveform');

            plot_spectrum(debug.rxSignal, cfg.sampleRate, 'Received spectrum');
        end

        if isfield(debug, 'sync') && isfield(debug.sync, 'correlationMetric')
            figure('Name', 'Preamble correlation metric');
            plot(debug.sync.correlationMetric);
            grid on;
            xlabel('Candidate packet start sample');
            ylabel('Normalized correlation');
            title('Preamble correlation metric');
        end

        if isfield(debug, 'rawDataSymbols')
            rawSymbols = debug.rawDataSymbols(:);
            rawSymbols = downsample_constellation(rawSymbols, cfg.maxConstellationPoints);
            figure('Name', 'Constellation before equalization');
            scatter(real(rawSymbols), imag(rawSymbols), 7, '.');
            axis equal;
            grid on;
            xlabel('In-phase');
            ylabel('Quadrature');
            title('Received constellation before equalization');
        end

        if isfield(debug, 'equalizedDataSymbols')
            eqSymbols = debug.equalizedDataSymbols(:);
            eqSymbols = downsample_constellation(eqSymbols, cfg.maxConstellationPoints);
            figure('Name', 'Constellation after equalization');
            scatter(real(eqSymbols), imag(eqSymbols), 7, '.');
            axis equal;
            grid on;
            xlabel('In-phase');
            ylabel('Quadrature');
            title('Received constellation after equalization');
        end

        if isfield(debug, 'pilotTrack')
            figure('Name', 'Continuous pilot amplitude');
            plot(debug.pilotTrack.amplitude.');
            grid on;
            xlabel('OFDM data symbol index');
            ylabel('Amplitude');
            title('Continuous pilots amplitude');
            legend(compose('active %d', cfg.contPilotActiveIndices), 'Location', 'best');

            figure('Name', 'Continuous pilot phase');
            plot(debug.pilotTrack.phase.');
            grid on;
            xlabel('OFDM data symbol index');
            ylabel('Phase (rad)');
            title('Continuous pilots phase');
            legend(compose('active %d', cfg.contPilotActiveIndices), 'Location', 'best');
        end
end
end

function plot_spectrum(x, sampleRate, plotTitle)
nFft = 2^nextpow2(min(numel(x), 2^18));
X = fftshift(fft(x(:), nFft));
fMHz = (-nFft/2:nFft/2-1).' * sampleRate / nFft / 1e6;
figure('Name', plotTitle);
plot(fMHz, 20 * log10(abs(X) / max(abs(X) + eps) + eps));
grid on;
xlabel('Frequency (MHz)');
ylabel('Magnitude (dB)');
title(plotTitle);
end

function symbols = downsample_constellation(symbols, maxPoints)
if numel(symbols) > maxPoints
    index = round(linspace(1, numel(symbols), maxPoints));
    symbols = symbols(index);
end
end
