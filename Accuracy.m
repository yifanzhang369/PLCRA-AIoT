function main()


    rng('shuffle');

   
    LArray = cell(1, 10);
    CArray = cell(1, 10);

    n = 1;

    for SNR = 0:1:4

        
        v = 3;                  
        fc = 9e8;                
        c = 3e8;                
        fd = v / c * fc;        
        Tc = 1 / fd;           
        Ts = Tc / 3;            
       

        symbolRate = 1e6;      
        num_trials = 1000;      

        
        bitsPerPhase = 200;

        if mod(bitsPerPhase, 2) ~= 0
            error('bitsPerPhase must be even for QPSK modulation.');
        end

        thresholds = 0:1:bitsPerPhase;

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        K1 = randi([0 1], 1, bitsPerPhase);   % 合法 PID key: K_i
        KR = randi([0 1], 1, bitsPerPhase);   % 合法 session key: K_{R,i}

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        L_distance = zeros(1, num_trials);

        for time = 1:num_trials
            gsmChan = stdchan('gsmEQx6', symbolRate, fd);

            L_distance(time) = ComputeLegitimateL1Distance( ...
                gsmChan, SNR, bitsPerPhase, K1, KR);
        end

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        C_distance = zeros(1, num_trials);

        for time = 1:num_trials
           
            gsmChan2 = stdchan('gsmEQx6', symbolRate, fd);  % RFS -> attacker
            gsmChan3 = stdchan('gsmEQx6', symbolRate, fd);   % attacker -> RFS

            C_distance(time) = ComputeAttackerL1Distance( ...
                gsmChan2, gsmChan3, SNR, bitsPerPhase, K1, KR);
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('n = 1');
        disp(n);
        disp(threshold_count_L_sum1);
        disp(threshold_count_C_sum1);

        LArray{n} = threshold_count_L_sum1;
        CArray{n} = threshold_count_C_sum1;

        n = n + 1;
    end

    save('SNR.mat', 'LArray', 'CArray');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function L1_distance = ComputeLegitimateL1Distance(gsmChan, SNR, bitsPerPhase, K1, KR)

    qpskMod = comm.QPSKModulator('BitInput', true);
    qpskDemod = comm.QPSKDemodulator('BitOutput', true);

    numSymbols = bitsPerPhase / 2;

    
    % Step 1: challenge
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    D1_bits = randi([0 1], 1, bitsPerPhase);
    D2_bits = randi([0 1], 1, bitsPerPhase);

    modSignal_D1 = qpskMod(D1_bits.');
    modSignal_D2 = qpskMod(D2_bits.');

    modSignal_K1 = qpskMod(K1.');
    modSignal_KR = qpskMod(KR.');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 2: RFS -> AIoT challenge transmission
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    H_RA_T1 = extractChannelResponse(gsmChan, numSymbols);
    H_RA_T2 = extractChannelResponse(gsmChan, numSymbols);

    y_A1 = modSignal_D1 .* H_RA_T1;
    y_A1 = awgn(y_A1, SNR, 'measured');

    y_A2 = modSignal_D2 .* H_RA_T2;
    y_A2 = awgn(y_A2, SNR, 'measured');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 3: AIoT backscatters two responses
    % b(T1)=K_i .* y_A(T1)
    % b(T2)=K_Ri .* y_A(T2)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    b1 = 0.25.*modSignal_K1 .* y_A1;
    b2 = 0.25.*modSignal_KR .* y_A2;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 4: AIoT -> RFS uplink backscatter
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    H_AR_T1 = extractChannelResponse(gsmChan, numSymbols);
    H_AR_T2 = extractChannelResponse(gsmChan, numSymbols);

    y_R1 = b1 .* H_AR_T1;
    y_R1 = awgn(y_R1, SNR, 'measured');

    y_R2 = b2 .* H_AR_T2;
    y_R2 = awgn(y_R2, SNR, 'measured');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 5: Ratio-based recovery at RFS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    z1 = safeDivide(y_R1, modSignal_D1);  
    z2 = safeDivide(y_R2, modSignal_D2);   

    ratioFeature = safeDivide(z1, z2);    
    Khat_symbol = ratioFeature .* modSignal_KR;

    
    receivedBits_Khat = qpskDemod(Khat_symbol);

    
    L1_distance = sum(abs(K1(:) - receivedBits_Khat(:)));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Impersonation attack
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function L1_distance = ComputeAttackerL1Distance(gsmChan2, gsmChan3, SNR, bitsPerPhase, K1, KR)

    qpskMod = comm.QPSKModulator('BitInput', true);
    qpskDemod = comm.QPSKDemodulator('BitOutput', true);

    numSymbols = bitsPerPhase / 2;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 1: RFS challenge 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    D1_bits = randi([0 1], 1, bitsPerPhase);
    D2_bits = randi([0 1], 1, bitsPerPhase);

    modSignal_D1 = qpskMod(D1_bits.');
    modSignal_D2 = qpskMod(D2_bits.');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 2
    % Random pseudo-PID key and random pseudo-session key
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    KE_bits  = randi([0 1], 1, bitsPerPhase);   % fake PID key
    KRE_bits = randi([0 1], 1, bitsPerPhase);   % fake session key

    modSignal_KE  = qpskMod(KE_bits.');
    modSignal_KRE = qpskMod(KRE_bits.');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 3: RFS -> attacker challenge transmission
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    H_RE_T1 = extractChannelResponse(gsmChan2, numSymbols);
    H_RE_T2 = extractChannelResponse(gsmChan2, numSymbols);

    y_E1 = modSignal_D1 .* H_RE_T1;
    y_E1 = awgn(y_E1, SNR, 'measured');

    y_E2 = modSignal_D2 .* H_RE_T2;
    y_E2 = awgn(y_E2, SNR, 'measured');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 4: attacker backscatters forged responses
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    forged_b1 = modSignal_KE  .* y_E1;
    forged_b2 = modSignal_KRE .* y_E2;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 5: attacker -> RFS uplink
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    H_ER_T1 = extractChannelResponse(gsmChan3, numSymbols);
    H_ER_T2 = extractChannelResponse(gsmChan3, numSymbols);

    y_R1 = forged_b1 .* H_ER_T1;
    y_R1 = awgn(y_R1, SNR, 'measured');

    y_R2 = forged_b2 .* H_ER_T2;
    y_R2 = awgn(y_R2, SNR, 'measured');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 6: RFS continues to perform recovery according to legitimate procedures.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    z1 = safeDivide(y_R1, modSignal_D1);
    z2 = safeDivide(y_R2, modSignal_D2);

    ratioFeature = safeDivide(z1, z2);

    modSignal_KR_legit = qpskMod(KR.');
    Khat_symbol = ratioFeature .* modSignal_KR_legit;

    receivedBits_Khat = qpskDemod(Khat_symbol);

    L1_distance = sum(abs(K1(:) - receivedBits_Khat(:)));
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