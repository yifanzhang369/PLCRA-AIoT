function main()
clc;
clear;
close all;

rng('shuffle');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SNR = 0;
num_trials = 1000;


LArray = cell(1, 3);
CArray = cell(1, 3);

n = 1;


v = 3;                  % 移动速度 (m/s)
fc = 9e8;               % 载波频率 (Hz)
c = 3e8;                % 光速
fd = v / c * fc;        % 多普勒频移
symbolRate = 1e6;       % 符号率

bit_length_xor = 20;
thresholds = 0:1:bit_length_xor;

L_distance = zeros(1, num_trials);
C_distance = zeros(1, num_trials);

for t = 1:num_trials
    [L_distance(t), C_distance(t)] = Baseline1_XOR_Replay_Current_EQ_Total(SNR, bit_length_xor);
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
    [L_distance(t), C_distance(t)] = Baseline2_CIR_Replay_Corr05(SNR, numTaps);
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
    [L_distance(t), C_distance(t)] = ThisWork_Replay_ChallengeExpired(SNR, bitsPerPhase, symbolRate, fd);
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


disp(threshold_count_L_sum1);
disp(threshold_count_C_sum1);

LArray{n} = threshold_count_L_sum1;
CArray{n} = threshold_count_C_sum1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 保存（保持原来的风格）
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('Replay.mat', 'LArray', 'CArray');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROC Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

legend('Baseline1', 'Baseline2', 'This work', ...
    'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight', 'bold');

hold off;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Baseline 1:
% challenge-response + XOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [L1_legit, L1_attack] = Baseline1_XOR_Replay_Current_EQ_Total(SNR, bit_length)

M = 4;
k = log2(M);

% 当前 challenge / secret key / response
challenge_bits = randi([0 1], 1, bit_length);
secret_key     = randi([0 1], 1, bit_length);
response_bits  = xor(challenge_bits, secret_key);

response_symbols = bi2de(reshape(response_bits, [], k));
response_mod     = pskmod(response_symbols, M, pi/M);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_LV = (randn + 1i*randn) / sqrt(2);   % legitimate -> verifier

y_legit = h_LV .* response_mod;
y_legit = awgn(y_legit, SNR, 'measured');

y_legit_eq = y_legit ./ (h_LV + 1e-12);

legit_demod = pskdemod(y_legit_eq, M, pi/M);
legit_bits  = de2bi(legit_demod, k);
legit_bits  = legit_bits(:)';

L1_legit = sum(abs(legit_bits - response_bits));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% legitimate -> attacker
h_LA = (randn + 1i*randn) / sqrt(2);
y_eve = h_LA .* response_mod;
y_eve = awgn(y_eve, SNR, 'measured');

% attacker amplify-and-forward
relayGain = 1;

% attacker -> verifier
h_AV = (randn + 1i*randn) / sqrt(2);
y_replay = h_AV .* (relayGain .* y_eve);
y_replay = awgn(y_replay, SNR, 'measured');

% verifier 用总等效信道 equalization
h_eq = h_AV * h_LA * relayGain;
y_replay_eq = y_replay ./ (h_eq + 1e-12);

replay_demod = pskdemod(y_replay_eq, M, pi/M);
replay_bits  = de2bi(replay_demod, k);
replay_bits  = replay_bits(:)';

L1_attack = sum(abs(replay_bits - response_bits));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Baseline 2:
% CIR-based authentication
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [L1_legit, L1_attack] = Baseline2_CIR_Replay_Corr05(SNR, numTaps)

rho_legit  = 0.98;
rho_attack = 0.8;
h_ref = (randn(numTaps,1) + 1i*randn(numTaps,1)) / sqrt(2);

h_legit = rho_legit * h_ref + sqrt(1 - rho_legit^2) * ...
          (randn(numTaps,1) + 1i*randn(numTaps,1)) / sqrt(2);

h_legit_est = awgn(h_legit, SNR, 'measured');

template_bits = CIR2Bits(h_ref);
legit_bits    = CIR2Bits(h_legit_est);

L1_legit = sum(abs(template_bits(:) - legit_bits(:)));


h_attack = rho_attack * h_ref + sqrt(1 - rho_attack^2) * ...
           (randn(numTaps,1) + 1i*randn(numTaps,1)) / sqrt(2);


h_attack_recorded = awgn(h_attack, SNR, 'measured');
h_attack_replayed = awgn(h_attack_recorded, SNR, 'measured');

attack_bits = CIR2Bits(h_attack_replayed);

L1_attack = sum(abs(template_bits(:) - attack_bits(:)));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLCRA-AIoT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [L1_legit, L1_attack] = ThisWork_Replay_ChallengeExpired(SNR, bitsPerPhase, symbolRate, fd)

qpskMod = comm.QPSKModulator('BitInput', true);
qpskDemod = comm.QPSKDemodulator('BitOutput', true);

numSymbols = bitsPerPhase / 2;
backscatterCoeff = 0.25;

% PID key + 当前共享 session key
K1 = randi([0 1], 1, bitsPerPhase);
KR = randi([0 1], 1, bitsPerPhase);

modK1 = qpskMod(K1.');
modKR = qpskMod(KR.');


D_cur_bits = randi([0 1], 1, bitsPerPhase);
D_old_bits = randi([0 1], 1, bitsPerPhase);

modD_cur = qpskMod(D_cur_bits.');
modD_old = qpskMod(D_old_bits.');

gsmChan_RA_cur = stdchan('gsmEQx6', symbolRate, fd);
H_RA_cur_1 = extractChannelResponse(gsmChan_RA_cur, numSymbols);
H_RA_cur_2 = extractChannelResponse(gsmChan_RA_cur, numSymbols);

yA_cur_1 = awgn(modD_cur .* H_RA_cur_1, SNR, 'measured');
yA_cur_2 = awgn(modD_cur .* H_RA_cur_2, SNR, 'measured');


b_cur_1 = backscatterCoeff .* modK1 .* yA_cur_1;
b_cur_2 = backscatterCoeff .* modKR .* yA_cur_2;


gsmChan_AR_cur = stdchan('gsmEQx6', symbolRate, fd);
H_AR_cur_1 = extractChannelResponse(gsmChan_AR_cur, numSymbols);
H_AR_cur_2 = extractChannelResponse(gsmChan_AR_cur, numSymbols);

yR_legit_1 = awgn(b_cur_1 .* H_AR_cur_1, SNR, 'measured');
yR_legit_2 = awgn(b_cur_2 .* H_AR_cur_2, SNR, 'measured');


z_legit_1 = safeDivide(yR_legit_1, modD_cur);
z_legit_2 = safeDivide(yR_legit_2, modD_cur);

ratio_legit = safeDivide(z_legit_1, z_legit_2);
Khat_legit_symbol = ratio_legit .* modKR;
Khat_legit_bits = qpskDemod(Khat_legit_symbol);

L1_legit = sum(abs(K1(:) - Khat_legit_bits(:)));


gsmChan_RA_old = stdchan('gsmEQx6', symbolRate, fd);
H_RA_old_1 = extractChannelResponse(gsmChan_RA_old, numSymbols);
H_RA_old_2 = extractChannelResponse(gsmChan_RA_old, numSymbols);

yA_old_1 = awgn(modD_old .* H_RA_old_1, SNR, 'measured');
yA_old_2 = awgn(modD_old .* H_RA_old_2, SNR, 'measured');


b_old_1 = backscatterCoeff .* modK1 .* yA_old_1;
b_old_2 = backscatterCoeff .* modKR .* yA_old_2;


gsmChan_AE_old = stdchan('gsmEQx6', symbolRate, fd);
H_AE_old_1 = extractChannelResponse(gsmChan_AE_old, numSymbols);
H_AE_old_2 = extractChannelResponse(gsmChan_AE_old, numSymbols);

yE_old_1 = awgn(b_old_1 .* H_AE_old_1, SNR, 'measured');
yE_old_2 = awgn(b_old_2 .* H_AE_old_2, SNR, 'measured');


gsmChan_ER = stdchan('gsmEQx6', symbolRate, fd);
H_ER_1 = extractChannelResponse(gsmChan_ER, numSymbols);
H_ER_2 = extractChannelResponse(gsmChan_ER, numSymbols);

yR_attack_1 = awgn(yE_old_1 .* H_ER_1, SNR, 'measured');
yR_attack_2 = awgn(yE_old_2 .* H_ER_2, SNR, 'measured');

% verifier 仍用当前 challenge D_cur 处理
z_attack_1 = safeDivide(yR_attack_1, modD_cur);
z_attack_2 = safeDivide(yR_attack_2, modD_cur);

ratio_attack = safeDivide(z_attack_1, z_attack_2);
Khat_attack_symbol = ratio_attack .* modKR;
Khat_attack_bits = qpskDemod(Khat_attack_symbol);

L1_attack = sum(abs(K1(:) - Khat_attack_bits(:)));

end


function bits = CIR2Bits(h)
bits = zeros(1, 2 * length(h));
for i = 1:length(h)
    bits(2*i-1) = real(h(i)) >= 0;
    bits(2*i)   = imag(h(i)) >= 0;
end
end

function complexGain = extractChannelResponse(chanObj, signalLength)

inputSignal = ones(signalLength, 1);
outputSignal = chanObj(inputSignal);

complexGain = zeros(signalLength, 1);
for i = 1:signalLength
    complexGain(i) = outputSignal(i) / inputSignal(i);
end

end


function out = safeDivide(a, b)
epsVal = 1e-12;
out = a ./ (b + epsVal);
end


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