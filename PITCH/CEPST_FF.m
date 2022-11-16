%% Calculating CEPSTRUM of signal  
function ff0 = CEPST_FF(y,fs)
%% PARAMETERS
freq_range = [80 350]; %frequency of pitch, needed for component at zero quefrency filtration 
quefency = fix([fs/freq_range(2) fs/freq_range(1)]);
%[y,a] = TSNR_FILT(y,fs,1,1);
overlap = 0.5;
% Window data
Window_length_sec = 30/1000; % in seconds

s_len = length(y);
window_length = fix(Window_length_sec*fs); % calcualting window length
window_type = hanning(window_length); %type of window function

% In some applications that process large amounts of data with fft, it is common to resize the input so that the number of samples is a power of 2.
FFT_length  = pow2(nextpow2(window_length));
 
%% CEPSTRUM
PC_OUT = zeros(length(y),1);
CC_OUT = zeros(length(y),1);
PHC_OUT = zeros(length(y),1);
rec = zeros(length(y),1);
% P_CEPST_L = zeros(FFT_length,1);
cnt = 1;
for p = 1:fix(window_length.*(1-overlap)):length(y)-FFT_length
    signal_fragment = y(p:p+window_length-1); % getiing part of signal for windowing
    sf_windowed = signal_fragment .* window_type; % applying winow function
    preemph = [1 -0.97]; sf_filt = filter(1,preemph,sf_windowed); % applying pre emphasis
    sfw_fft = fft(sf_filt, FFT_length); % calculate fft of windowed part
    sq_phase = angle(sfw_fft); % Calculate phase of signal
    %% NORMAL CEPSTRUM
    N_CEPST = ifft(log(abs(sfw_fft)),FFT_length);
    %% POWER CEPSTRUM
    P_CEPST = abs(ifft(log(abs(sfw_fft).^2),FFT_length)).^2;
    PC_OUT(p:p+FFT_length-1) = PC_OUT(p:p+FFT_length-1) + P_CEPST*overlap;
    %% COMPLEX CEPSTRUM
    C_CEPST = ifft((log(abs(sfw_fft))+1i*sq_phase),FFT_length);
    CC_OUT(p:p+FFT_length-1) = CC_OUT(p:p+FFT_length-1) + C_CEPST*overlap; 
    recc = ifft(exp(fft(C_CEPST)));
    rec(p:p+FFT_length-1) = rec(p:p+FFT_length-1) + recc*overlap;
    %% PHASE CEPSTRUM
    PH_CEPST = (C_CEPST-recc).^2;
    PHC_OUT(p:p+FFT_length-1) = PHC_OUT(p:p+FFT_length-1) + PH_CEPST/(1/overlap);
    %% LOW TIME LIFTERING 
    N = 150;
    ltl = [ones(N,1); zeros(length(N_CEPST)-N,1)];
    P_CEPST_L = ifft(exp(abs(fft(ltl .* N_CEPST, FFT_length))));
    %[pks,lock] = findpeaks(P_CEPST_L(yy-60:yy));
    %locks{cnt,:} = lock+yy-60;
    %% PITCH ESTIMATION
    [xx,yy] = max(P_CEPST(quefency(1):quefency(2)));
    E(cnt) = xx/P_CEPST(1);
    if xx/P_CEPST(1)>= 0.0002
        ff0(cnt) = fs/(yy+quefency(1));
    else
        ff0(cnt) = 0;
    end
    cnt = cnt+1;
end
ff0 = ff0(ff0~=0);
end