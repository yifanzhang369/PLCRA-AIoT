function main()
clc;
clear;
close all;

rng('shuffle');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SNR sweep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SNR_dB = 0:5:25;

% Leak information curves
I_baseline1 = zeros(1, length(SNR_dB));   % Challenge-response + XOR
I_baseline2 = zeros(1, length(SNR_dB));   % CIR-based authentication
I_thiswork  = zeros(1, length(SNR_dB));   % Proposed ratio-based scheme

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_trials = 1000;          % Monte Carlo trials per SNR
Fc = 900e6;                % Carrier frequency
symbolRate = 1e6;          % Symbol rate for time-varying channel
v = 3;                     % Mobility speed (m/s)
c = 3e8;                   % Speed of light
fd = v / c * Fc;           % Doppler frequency

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Baseline1 parameters (XOR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bit_length_xor = 20;
M = 4;                     % QPSK
k = log2(M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Baseline2 parameters (CIR-based authentication)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numTaps = 10;              % 10 complex taps -> 20 bits after sign quantization
eve_distance = 0.5;        % distance between Eve and legit device for spatial correlation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This work parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bitsPerPhase = 200;        % PID key length
backscatterCoeff = 0.25;   % reflection coefficient

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(SNR_dB)
    SNR = SNR_dB(i);

    % Baseline1: challenge-response + XOR
    I_baseline1(i) = Baseline1_XOR_LeakInfo(SNR, num_trials, bit_length_xor, M, k);

    % Baseline2: CIR-based authentication
    I_baseline2(i) = Baseline2_CIR_LeakInfo(SNR, num_trials, numTaps, Fc, eve_distance);

    % This work: ratio-based two-response backscatter authentication
    I_thiswork(i) = ThisWork_LeakInfo(SNR, num_trials, bitsPerPhase, ...
                                      symbolRate, fd, backscatterCoeff);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
markers = {'o', 's', 'd', 'v', 'p', '*'};
[~, all_colors] = GetColors();

figure;
hold on;

plot(SNR_dB, I_baseline1, 'Color', all_colors(5, :), ...
    'Marker', markers{1}, 'LineStyle', '-', 'LineWidth', 2.5);

plot(SNR_dB, I_baseline2, 'Color', all_colors(4, :), ...
    'Marker', markers{2}, 'LineStyle', '-', 'LineWidth', 2.5);

plot(SNR_dB, I_thiswork, 'Color', all_colors(3, :), ...
    'Marker', markers{3}, 'LineStyle', '-', 'LineWidth', 2.5);

xlabel('SNR (dB)', 'FontWeight', 'bold', 'FontSize', 20, 'FontName', 'Times New Roman');
ylabel('Leak Information (per bit)', 'FontWeight', 'bold', 'FontSize', 20, 'FontName', 'Times New Roman');

xticks(0:5:25);
yticks(0:0.2:1);

ax = gca;
ax.XMinorTick = 'off';
ax.YMinorTick = 'off';
xlim([0 25]);
ylim([0 1]);
ax.Box = 'on';
ax.XColor = 'k';
ax.YColor = 'k';

ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
ax.XAxis.FontName = 'Times New Roman';
ax.YAxis.FontName = 'Times New Roman';

grid on;

legend('Trad', 'STA', 'This work', ...
    'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight', 'bold');

hold off;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Baseline1:
% Challenge-response + XOR
%
% Eavesdropper intercepts both the plaintext challenge and the XOR-masked
% response, then guesses the secret key by XORing them.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function I_avg = Baseline1_XOR_LeakInfo(SNR, num_trials, bit_length, M, k)

I_acc = 0;

for t = 1:num_trials
    % Random challenge and secret key
    challenge_bits = randi([0 1], 1, bit_length);
    secret_key = randi([0 1], 1, bit_length);

    % Legitimate response: challenge XOR secret_key
    response_bits = xor(challenge_bits, secret_key);

    % QPSK modulation
    challenge_symbols = bi2de(reshape(challenge_bits, [], k));
    challenge_modulated = pskmod(challenge_symbols, M, pi/M);

    response_symbols = bi2de(reshape(response_bits, [], k));
    response_modulated = pskmod(response_symbols, M, pi/M);

    % Eve intercepts challenge and response over AWGN
    challenge_intercepted = awgn(challenge_modulated, SNR, 'measured');
    response_intercepted  = awgn(response_modulated,  SNR, 'measured');

    % Demodulation at Eve
    challenge_demodulated = pskdemod(challenge_intercepted, M, pi/M);
    response_demodulated  = pskdemod(response_intercepted,  M, pi/M);

    challenge_eavesdropped = de2bi(challenge_demodulated, k);
    challenge_eavesdropped = challenge_eavesdropped(:)';

    response_eavesdropped = de2bi(response_demodulated, k);
    response_eavesdropped = response_eavesdropped(:)';

    % Eve guesses the key by XOR
    guessed_key = xor(challenge_eavesdropped, response_eavesdropped);

    % Leak information
    I_acc = I_acc + mutual_information2(guessed_key, secret_key, 2);
end

I_avg = I_acc / num_trials;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Baseline2:
% CIR-based authentication
%
% Legitimate authentication feature is derived from the CIR.
% Eve tries to infer the legitimate CIR-based feature from its own observed
% correlated CIR. Leakage is measured by MI between legitimate CIR bits and
% Eve's CIR bits.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function I_avg = Baseline2_CIR_LeakInfo(SNR, num_trials, numTaps, Fc, eve_distance)

c = 3e8;
lambda = c / Fc;

% Simple spatial correlation model
rho = abs(besselj(0, 2*pi*eve_distance/lambda));
rho = min(max(rho, 0), 0.999);

I_acc = 0;

for t = 1:num_trials
    % Legitimate CIR
    h_legit = (randn(numTaps,1) + 1i*randn(numTaps,1)) / sqrt(2);

    % Eve CIR correlated with legitimate CIR
    h_eve = rho * h_legit + sqrt(1 - rho^2) * ...
            (randn(numTaps,1) + 1i*randn(numTaps,1)) / sqrt(2);

    % Noisy estimation at both sides
    h_legit_est = awgn(h_legit, SNR, 'measured');
    h_eve_est   = awgn(h_eve,   SNR, 'measured');

    % Quantize CIR to authentication bits
    legit_bits = CIR2Bits(h_legit_est);
    eve_bits   = CIR2Bits(h_eve_est);

    % Leak information
    I_acc = I_acc + mutual_information2(eve_bits, legit_bits, 2);
end

I_avg = I_acc / num_trials;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This work:
% Single challenge + two backscatter responses + ratio-based masking
%
% Eve intercepts the two challenge blocks and the two backscatter responses.
% After removing the challenge, Eve can only obtain a masked ratio feature
% that depends on K_i / K_Ri, rather than K_i itself.
% Leakage is measured between Eve's extracted observation bits and the true
% PID key bits.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function I_avg = ThisWork_LeakInfo(SNR, num_trials, bitsPerPhase, symbolRate, fd, backscatterCoeff)

qpskMod = comm.QPSKModulator('BitInput', true);
qpskDemod = comm.QPSKDemodulator('BitOutput', true);

numSymbols = bitsPerPhase / 2;
I_acc = 0;

for t = 1:num_trials
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fresh session keys
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    K1_bits = randi([0 1], 1, bitsPerPhase);   % PID key
    KR_bits = randi([0 1], 1, bitsPerPhase);   % session key

    D1_bits = randi([0 1], 1, bitsPerPhase);   % challenge block 1
    D2_bits = randi([0 1], 1, bitsPerPhase);   % challenge block 2

    modK1 = qpskMod(K1_bits.');
    modKR = qpskMod(KR_bits.');
    modD1 = qpskMod(D1_bits.');
    modD2 = qpskMod(D2_bits.');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % RFS -> AIoT challenge transmission
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gsmChan_RA = stdchan('gsmEQx6', symbolRate, fd);
    H_RA_1 = extractChannelResponse(gsmChan_RA, numSymbols);
    H_RA_2 = extractChannelResponse(gsmChan_RA, numSymbols);

    yA1 = awgn(modD1 .* H_RA_1, SNR, 'measured');
    yA2 = awgn(modD2 .* H_RA_2, SNR, 'measured');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % AIoT generates two backscatter responses
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    b1 = backscatterCoeff .* modK1 .* yA1;   % response carrying PID key
    b2 = backscatterCoeff .* modKR .* yA2;   % response carrying session key

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % AIoT -> Eve observation channels
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gsmChan_AE = stdchan('gsmEQx6', symbolRate, fd);
    H_AE_1 = extractChannelResponse(gsmChan_AE, numSymbols);
    H_AE_2 = extractChannelResponse(gsmChan_AE, numSymbols);

    yE1 = awgn(b1 .* H_AE_1, SNR, 'measured');
    yE2 = awgn(b2 .* H_AE_2, SNR, 'measured');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Eve knows/decodes D1 and D2, so it removes challenge symbols
    % Then it forms a ratio feature.
    % This reveals only a masked feature related to K_i / K_Ri.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    zE1 = safeDivide(yE1, modD1);
    zE2 = safeDivide(yE2, modD2);

    ratioFeatureE = safeDivide(zE1, zE2);

    % Eve demodulates the extracted observation
    observed_bits = qpskDemod(ratioFeatureE);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Leak information between Eve's observation and the true PID key
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    I_acc = I_acc + mutual_information2(observed_bits, K1_bits, 2);
end

I_avg = I_acc / num_trials;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert CIR to bits:
% Each complex tap gives two bits:
%   bit 1 = sign(real(h))
%   bit 2 = sign(imag(h))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bits = CIR2Bits(h)

bits = zeros(1, 2*length(h));

for i = 1:length(h)
    bits(2*i-1) = real(h(i)) >= 0;
    bits(2*i)   = imag(h(i)) >= 0;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract time-varying channel response from stdchan object
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function complexGain = extractChannelResponse(chanObj, signalLength)

inputSignal = ones(signalLength, 1);
outputSignal = chanObj(inputSignal);

complexGain = zeros(signalLength, 1);
for i = 1:signalLength
    complexGain(i) = outputSignal(i) / inputSignal(i);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Safe divide
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = safeDivide(a, b)
epsVal = 1e-12;
out = a ./ (b + epsVal);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mutual information (normalized)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function NMI = mutual_information2(x, y, num_bins)
    % x, y: two input vectors
    % num_bins: number of bins for discretization
    % output NMI: normalized mutual information in [0,1]

    if length(x) ~= length(y)
        error('Input vectors x and y must have the same length.');
    end

    if nargin < 3
        num_bins = ceil(log2(length(x)) + 1);
    end

    x = x(:);
    y = y(:);

    joint_hist = histcounts2(x, y, num_bins, 'Normalization', 'probability');

    px = sum(joint_hist, 2);
    py = sum(joint_hist, 1);

    MI = 0;
    for i = 1:num_bins
        for j = 1:num_bins
            if joint_hist(i, j) > 0 && px(i) > 0 && py(j) > 0
                MI = MI + joint_hist(i, j) * log2(joint_hist(i, j) / (px(i) * py(j)));
            end
        end
    end

    Hx = -sum(px(px > 0) .* log2(px(px > 0)));
    Hy = -sum(py(py > 0) .* log2(py(py > 0)));

    if Hx == 0 || Hy == 0
        NMI = 0;
        return;
    end

    NMI = MI / max(Hx, Hy);

    % numerical protection
    NMI = max(0, min(1, real(NMI)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Color utilities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [all_themes, all_colors] = GetColors()
c_string{1} = {'#FD6D5A', '#FEB40B', '#6DC354', '#994487', '#518CD8', '#443295'};
c_string{2} = {'#264653', '#2A9D8F', '#E9C46A', '#F4A261', '#E76F51', '#253777'};
c_string{3} = {'#C1C976', '#C8A9A1', '#FEC2E4', '#77CCE0', '#FFD372', '#F88078'};
c_string{4} = {'#104FFF', '#2FD151', '#64C7B8', '#FF1038', '#45CAFF', '#B913FF'};
c_string{5} = {'#4C87D6', '#F38562', '#F2B825', '#D4C114', '#88B421', '#199FE0'};
c_string{6} = {'#037CD2', '#00AAAA', '#927FD3', '#E54E5D', '#EAA700', '#F57F4B'};

m_count = length(c_string);
c_count = sum(cellfun(@length, c_string));
all_colors = nan(c_count, 3);
all_themes = cell(m_count, 1);
idx = 1;

for i = 1:length(c_string)
    color = nan(length(c_string{i}), 3);
    for j = 1:length(c_string{i})
        rgb = Hex2RGB(c_string{i}{j});
        all_colors(idx, :) = RGB2MatlabColor(rgb);
        color(j, :) = all_colors(idx, :);
        idx = idx + 1;
    end
    all_themes{i} = color;
end
end

function c = Hex2RGB(str)
clist = '0123456789ABCDEF';
nums = zeros(1, 6);
for i = 1:6
    nums(i) = find(str(i + 1) == clist);
end
c = zeros(1, 3);
for i = 1:3
    c(i) = 16 * (nums(2*i-1) - 1) + (nums(2*i) - 1);
end
end

function c = RGB2MatlabColor(rgb)
c = round(100 * rgb / 255) / 100;
end