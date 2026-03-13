function main()
clc;
clear;
close all;

rng('shuffle');


SNR = 0;
num_trials = 1000;


LArray = cell(1, 3);
CArray = cell(1, 3);

n = 1;


v = 3;                 
fc = 9e8;               
c = 3e8;                
fd = v / c * fc;        
symbolRate = 1e6;    

bit_length_xor = 20;
thresholds = 0:1:bit_length_xor;

L_distance = zeros(1, num_trials);
C_distance = zeros(1, num_trials);

for t = 1:num_trials
    [L_distance(t), C_distance(t)] = Baseline1_XOR_Impersonation_GSM( ...
        SNR, bit_length_xor, symbolRate, fd);
end

threshold_count_L_sum1 = zeros(1, length(thresholds));
threshold_count_C_sum1 = zeros(1, length(thresholds));

for j = 1:length(thresholds)
    delta = thresholds(j);
    for i = 1:num_trials
        if L_distance(i) <= delta
            threshold_count_L_sum1(j) = threshold_count_L_sum1(j) + 1;
        end
        if C_distance(i) <= delta
            threshold_count_C_sum1(j) = threshold_count_C_sum1(j) + 1;
        end
    end
end

disp('n = 1');
disp(n);
disp(threshold_count_L_sum1);
disp(threshold_count_C_sum1);

LArray{n} = threshold_count_L_sum1;
CArray{n} = threshold_count_C_sum1;
n = n + 1;

numTaps = 20;
cir_bit_length = 2 * numTaps;
thresholds = 0:1:cir_bit_length;

L_distance = zeros(1, num_trials);
C_distance = zeros(1, num_trials);

for t = 1:num_trials
    [L_distance(t), C_distance(t)] = Baseline2_CIR_Impersonation_GSM( ...
        SNR, numTaps, symbolRate, fd);
end

threshold_count_L_sum1 = zeros(1, length(thresholds));
threshold_count_C_sum1 = zeros(1, length(thresholds));

for j = 1:length(thresholds)
    delta = thresholds(j);
    for i = 1:num_trials
        if L_distance(i) <= delta
            threshold_count_L_sum1(j) = threshold_count_L_sum1(j) + 1;
        end
        if C_distance(i) <= delta
            threshold_count_C_sum1(j) = threshold_count_C_sum1(j) + 1;
        end
    end
end

disp('n = 1');
disp(n);
disp(threshold_count_L_sum1);
disp(threshold_count_C_sum1);

LArray{n} = threshold_count_L_sum1;
CArray{n} = threshold_count_C_sum1;
n = n + 1;


bitsPerPhase = 200;
thresholds = 0:1:bitsPerPhase;

L_distance = zeros(1, num_trials);
C_distance = zeros(1, num_trials);

for t = 1:num_trials
    [L_distance(t), C_distance(t)] = ThisWork_Impersonation( ...
        SNR, bitsPerPhase, symbolRate, fd);
end

threshold_count_L_sum1 = zeros(1, length(thresholds));
threshold_count_C_sum1 = zeros(1, length(thresholds));

for j = 1:length(thresholds)
    delta = thresholds(j);
    for i = 1:num_trials
        if L_distance(i) <= delta
            threshold_count_L_sum1(j) = threshold_count_L_sum1(j) + 1;
        end
        if C_distance(i) <= delta
            threshold_count_C_sum1(j) = threshold_count_C_sum1(j) + 1;
        end
    end
end

disp('n = 1');
disp(n);
disp(threshold_count_L_sum1);
disp(threshold_count_C_sum1);

LArray{n} = threshold_count_L_sum1;
CArray{n} = threshold_count_C_sum1;
n = n + 1;

markers = {'o', 's', 'd', 'v', 'p', '*'};
[~, all_colors] = GetColors();

figure;
hold on;

plot(CArray{1}/num_trials, LArray{1}/num_trials, ...
    'Color', all_colors(5, :), 'Marker', markers{1}, ...
    'LineStyle', '-', 'LineWidth', 2.5);

plot(CArray{2}/num_trials, LArray{2}/num_trials, ...
    'Color', all_colors(4, :), 'Marker', markers{2}, ...
    'LineStyle', '-', 'LineWidth', 2.5);

plot(CArray{3}/num_trials, LArray{3}/num_trials, ...
    'Color', all_colors(3, :), 'Marker', markers{3}, ...
    'LineStyle', '-', 'LineWidth', 2.5);

plot([0 1], [0 1], '--k', 'LineWidth', 1.5);

xlabel('False Acceptance Rate', 'FontWeight', 'bold', ...
    'FontSize', 20, 'FontName', 'Times New Roman');
ylabel('Authentication Rate', 'FontWeight', 'bold', ...
    'FontSize', 20, 'FontName', 'Times New Roman');

xticks(0:0.2:1);
yticks(0:0.2:1);

ax = gca;
ax.XMinorTick = 'off';
ax.YMinorTick = 'off';
xlim([0 1]);
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

legend('Tradition CR', 'Related PLA', 'PLCRA-AIoT', ...
    'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight', 'bold');

hold off;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Baseline 1
% challenge-response + XOR
%
% Counterfeiting attack via eavesdropping:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [L1_legit, L1_attack] = Baseline1_XOR_Impersonation_GSM(SNR, bit_length, symbolRate, fd)

M = 4;
k = log2(M);


SNR_EVE = SNR + 10;

secret_key = randi([0 1], 1, bit_length);


challenge_prev = randi([0 1], 1, bit_length);
response_prev  = xor(challenge_prev, secret_key);

challenge_prev_symbols = bi2de(reshape(challenge_prev, [], k));
challenge_prev_mod     = pskmod(challenge_prev_symbols, M, pi/M);

response_prev_symbols = bi2de(reshape(response_prev, [], k));
response_prev_mod     = pskmod(response_prev_symbols, M, pi/M);

challenge_prev_hat = GSM_QPSK_Receive(challenge_prev_mod, SNR_EVE, symbolRate, fd, k);
response_prev_hat  = GSM_QPSK_Receive(response_prev_mod,  SNR_EVE, symbolRate, fd, k);

challenge_prev_hat = challenge_prev_hat(1:bit_length);
response_prev_hat  = response_prev_hat(1:bit_length);

key_est_at_eve = xor(challenge_prev_hat, response_prev_hat);


challenge_cur = randi([0 1], 1, bit_length);
response_cur  = xor(challenge_cur, secret_key);

challenge_cur_symbols = bi2de(reshape(challenge_cur, [], k));
challenge_cur_mod     = pskmod(challenge_cur_symbols, M, pi/M);

response_cur_symbols = bi2de(reshape(response_cur, [], k));
response_cur_mod     = pskmod(response_cur_symbols, M, pi/M);


legit_bits = GSM_QPSK_Receive(response_cur_mod, SNR, symbolRate, fd, k);
legit_bits = legit_bits(1:bit_length);


L1_legit = sum(abs(legit_bits - response_cur));

challenge_cur_hat = GSM_QPSK_Receive(challenge_cur_mod, SNR_EVE, symbolRate, fd, k);
challenge_cur_hat = challenge_cur_hat(1:bit_length);

forged_response = xor(challenge_cur_hat, key_est_at_eve);

forged_symbols = bi2de(reshape(forged_response, [], k));
forged_mod     = pskmod(forged_symbols, M, pi/M);

attack_bits = GSM_QPSK_Receive(forged_mod, SNR, symbolRate, fd, k);
attack_bits = attack_bits(1:bit_length);

L1_attack = sum(abs(attack_bits - response_cur));

end

function [L1_legit, L1_attack] = Baseline2_CIR_Impersonation_GSM(SNR, numTaps, symbolRate, fd)

rho_legit = 0.98;

gsmChan_ref = stdchan('gsmEQx6', symbolRate, fd);
H_ref_est = extractChannelResponse(gsmChan_ref, numTaps);
H_ref_sig = extractChannelResponse(gsmChan_ref, numTaps);

H_ref_feature = 0.5 * (H_ref_est + H_ref_sig);
template_bits = CIR2Bits(H_ref_feature);

gsmChan_legit = stdchan('gsmEQx6', symbolRate, fd);
H_legit_est_base = extractChannelResponse(gsmChan_legit, numTaps);
H_legit_sig_base = extractChannelResponse(gsmChan_legit, numTaps);

H_legit_est = rho_legit * H_ref_est + sqrt(1 - rho_legit^2) * H_legit_est_base;
H_legit_sig = rho_legit * H_ref_sig + sqrt(1 - rho_legit^2) * H_legit_sig_base;

H_legit_sig = awgn(H_legit_sig, SNR, 'measured');

H_legit_feature = 0.5 * (H_legit_est + H_legit_sig);
legit_bits = CIR2Bits(H_legit_feature);

L1_legit = sum(abs(template_bits(:) - legit_bits(:)));

gsmChan_attack = stdchan('gsmEQx6', symbolRate, fd);
H_attack_est_base = extractChannelResponse(gsmChan_attack, numTaps);
H_attack_sig_base = extractChannelResponse(gsmChan_attack, numTaps);

theta = 2*pi*rand;
rotationFactor = exp(1j * theta);

H_attack_est = rotationFactor .* H_attack_est_base;
H_attack_sig = rotationFactor .* H_attack_sig_base;
H_attack_sig = awgn(H_attack_sig, SNR, 'measured');

H_attack_feature = 0.5 * (H_attack_est + H_attack_sig);
attack_bits = CIR2Bits(H_attack_feature);

L1_attack = sum(abs(template_bits(:) - attack_bits(:)));

end

function [L1_legit, L1_attack] = ThisWork_Impersonation(SNR, bitsPerPhase, symbolRate, fd)

qpskMod = comm.QPSKModulator('BitInput', true);
qpskDemod = comm.QPSKDemodulator('BitOutput', true);

numSymbols = bitsPerPhase / 2;
backscatterCoeff = 0.25;

K1 = randi([0 1], 1, bitsPerPhase);
KR = randi([0 1], 1, bitsPerPhase);

modK1 = qpskMod(K1.');
modKR = qpskMod(KR.');

D_bits = randi([0 1], 1, bitsPerPhase);
modD   = qpskMod(D_bits.');

gsmChan_RA = stdchan('gsmEQx6', symbolRate, fd);
H_RA_1 = extractChannelResponse(gsmChan_RA, numSymbols);
H_RA_2 = extractChannelResponse(gsmChan_RA, numSymbols);

yA_1 = awgn(modD .* H_RA_1, SNR, 'measured');
yA_2 = awgn(modD .* H_RA_2, SNR, 'measured');

b1 = backscatterCoeff .* modK1 .* yA_1;
b2 = backscatterCoeff .* modKR .* yA_2;

gsmChan_AR = stdchan('gsmEQx6', symbolRate, fd);
H_AR_1 = extractChannelResponse(gsmChan_AR, numSymbols);
H_AR_2 = extractChannelResponse(gsmChan_AR, numSymbols);

yR_legit_1 = awgn(b1 .* H_AR_1, SNR, 'measured');
yR_legit_2 = awgn(b2 .* H_AR_2, SNR, 'measured');

z_legit_1 = safeDivide(yR_legit_1, modD);
z_legit_2 = safeDivide(yR_legit_2, modD);

ratio_legit = safeDivide(z_legit_1, z_legit_2);
Khat_legit_symbol = ratio_legit .* modKR;
Khat_legit_bits = qpskDemod(Khat_legit_symbol);

L1_legit = sum(abs(K1(:) - Khat_legit_bits(:)));

KE  = randi([0 1], 1, bitsPerPhase);   % fake PID key
KRE = randi([0 1], 1, bitsPerPhase);   % fake session key

modKE  = qpskMod(KE.');
modKRE = qpskMod(KRE.');

gsmChan_RE = stdchan('gsmEQx6', symbolRate, fd);
H_RE_1 = extractChannelResponse(gsmChan_RE, numSymbols);
H_RE_2 = extractChannelResponse(gsmChan_RE, numSymbols);

yE_1 = awgn(modD .* H_RE_1, SNR, 'measured');
yE_2 = awgn(modD .* H_RE_2, SNR, 'measured');

forged_b1 = backscatterCoeff .* modKE  .* yE_1;
forged_b2 = backscatterCoeff .* modKRE .* yE_2;

gsmChan_ER = stdchan('gsmEQx6', symbolRate, fd);
H_ER_1 = extractChannelResponse(gsmChan_ER, numSymbols);
H_ER_2 = extractChannelResponse(gsmChan_ER, numSymbols);

yR_attack_1 = awgn(forged_b1 .* H_ER_1, SNR, 'measured');
yR_attack_2 = awgn(forged_b2 .* H_ER_2, SNR, 'measured');

z_attack_1 = safeDivide(yR_attack_1, modD);
z_attack_2 = safeDivide(yR_attack_2, modD);

ratio_attack = safeDivide(z_attack_1, z_attack_2);
Khat_attack_symbol = ratio_attack .* modKR;
Khat_attack_bits = qpskDemod(Khat_attack_symbol);

L1_attack = sum(abs(K1(:) - Khat_attack_bits(:)));

end

function rx_bits = GSM_QPSK_Receive(tx_mod, SNR, symbolRate, fd, k)

numSymbols = length(tx_mod);

gsmChan = stdchan('gsmEQx6', symbolRate, fd);

% estimation phase
H_est = extractChannelResponse(gsmChan, numSymbols);

% authentication/data phase
H_sig = extractChannelResponse(gsmChan, numSymbols);

% receive
y = awgn(tx_mod .* H_sig, SNR, 'measured');

% equalization by estimated channel
y_eq = safeDivide(y, H_est);

% demodulation
rx_demod = pskdemod(y_eq, 4, pi/4);
rx_bits = de2bi(rx_demod, k);
rx_bits = rx_bits(:)';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CIR -> bits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bits = CIR2Bits(h)
bits = zeros(1, 2 * length(h));
for i = 1:length(h)
    bits(2*i-1) = real(h(i)) >= 0;
    bits(2*i)   = imag(h(i)) >= 0;
end
end


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
function out = safeDivide(a, b)
epsVal = 1e-12;
out = a ./ (b + epsVal);
end

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