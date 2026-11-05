clear; clc; close all;

% Parameters
GroupNumber = 16;
alpha = GroupNumber/18;
mu = 0.00022676;
Fs = 44100; % sample
fc = 4000; % carrier
Tsym = 2.2676e-3;
Ns = round(Tsym*Fs);
Ts = 1/Fs;
mu_samp = round(mu*Fs);

M = 4;
bitsPerSym = log2(M);
EbN0dB = 0:9;
Nbits = 1e6;
Nsym = Nbits/bitsPerSym;

% Pulse generation
n = 0:Ns-1;
p = sin(pi*n/Ns);
p = p / sqrt(sum(p.^2));
mf = fliplr(p);

bitsTx = randi([0 1], Nbits, 1);
symTx = map_qpsk(bitsTx);
I = real(symTx);
Q = imag(symTx);

% Baseband signal generation
Iup = zeros(Nsym*Ns,1);
Qup = zeros(Nsym*Ns,1);
Iup(1:Ns:end) = I;
Qup(1:Ns:end) = Q;
Ibb = conv(Iup, p);
Qbb = conv(Qup, p);

% Carrier generation and modulation
ntx = 0:length(Ibb)-1;
carrierCos = cos(2*pi*fc*ntx/Fs).';
carrierSin = sin(2*pi*fc*ntx/Fs).';

% s(t) = I*cos - Q*sin
sTx = Ibb .* carrierCos - Qbb .* carrierSin;

% Channel
h = [sqrt(1-alpha^2), zeros(1, mu_samp-1), alpha];
rClean = conv(sTx, h);

% Eb calculation
Etx_total = sum(sTx.^2) / Fs;
Eb = Etx_total / Nbits;

sampleOffset = Ns;

BER_sim = zeros(size(EbN0dB));
BER_theory = q_function(sqrt(2*10.^(EbN0dB/10)));

nrx = 0:length(rClean)-1;
loCos = 2*cos(2*pi*fc*nrx/Fs).';
loSin = -2*sin(2*pi*fc*nrx/Fs).';

for k = 1:length(EbN0dB)
    EbN0dB_linear = 10^(EbN0dB(k)/10);
    N0 = Eb / EbN0dB_linear;
    sigma2 = N0*Fs/2;
    noise = sqrt(sigma2) * randn(size(rClean));
    r = rClean + noise;

    % Downconversion
    rI = r .* loCos;
    rQ = r .* loSin;

    % Matched filter
    zI = conv(rI, mf);
    zQ = conv(rQ, mf);

    sampleIdx = sampleOffset:Ns:(sampleOffset + (Nsym-1)*Ns);
    z = zI(sampleIdx) + 1j*zQ(sampleIdx);

    bits_receive = demap_qpsk(z);
    BER_sim(k) = mean(bits_receive ~= bitsTx);

    fprintf('Eb/N0 = %2d dB, BER = %.6e\n', EbN0dB(k), BER_sim(k));
end

% Plot
figure;
semilogy(EbN0dB, BER_sim); hold on;
semilogy(EbN0dB, BER_theory);
grid on;
xlabel('E_b/N_0 (dB)', 'FontSize', 18);
ylabel('BER', 'FontSize', 18);
title('Task 1 BER', 'FontSize', 18);
legend('Simulation', 'Theory', 'FontSize', 18);
ylim([1e-5 1e-1]);
set(gca, 'FontSize', 18)

figure;
plot(zI(sampleIdx), zQ(sampleIdx),'.');
title('Constellation', 'FontSize', 18);
set(gca, 'FontSize', 18)

% Functions
function s = map_qpsk(bits)
    bits = bits(:);
    b = reshape(bits, 2, []).';
    s = zeros(size(b,1),1);

    idx = (b(:,1)==0 & b(:,2)==0); s(idx) = (1 + 1j)/sqrt(2);  % 00
    idx = (b(:,1)==0 & b(:,2)==1); s(idx) = (-1 + 1j)/sqrt(2);  % 01
    idx = (b(:,1)==1 & b(:,2)==1); s(idx) = (-1 - 1j)/sqrt(2);  % 11
    idx = (b(:,1)==1 & b(:,2)==0); s(idx) = (1 - 1j)/sqrt(2);  % 10
end

function bits = demap_qpsk(z)
    z = z(:);
    bits = zeros(2*length(z),1);
    for m = 1:length(z)
        if real(z(m)) >= 0 && imag(z(m)) >= 0
            pair = [0;0];
        elseif real(z(m)) < 0 && imag(z(m)) >= 0
            pair = [0;1];
        elseif real(z(m)) < 0 && imag(z(m)) < 0
            pair = [1;1];
        else
            pair = [1;0];
        end
        bits(2*m-1:2*m) = pair;
    end
end

function y = q_function(x)
    y = 0.5*erfc(x./sqrt(2));
end
