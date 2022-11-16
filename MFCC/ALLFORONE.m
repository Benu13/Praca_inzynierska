function [f0,c_coef,d_coef, formants, hhb] = ALLFORONE(y,fs,Mel_coefs,d_on,yoke)
%ALLFORONE Summary of this function goes here
%   Detailed explanation goes here
%% MFCC, Mel-frequency cepstral coefficients
% Najpierw obliczymy ton podstawowy z wykorzystaniem ACF, przy okazji
% wykorzystamy decyzję czy sygnał jest voiced/unvoiced żeby stworzyć
% kombinację samych ramek które zostały określone jako voiced co pozwoli
% nam pominąć problem szumu. Na poswstałym zlepku będziemy badać MFCC. 
step = 1;
f0 = [];
CEPST_COEFS = [];
FORMANTS_C = [];
hhb = [];
%% TSNR_filtration
%[~,y] = TSNR_FILT(y,fs);

B = [1 -0.98];
%% Pitch
while (isempty(f0) && step <= 4)
    ys = y;
    % Get pitch
    [f0, ys] = ACF_FF(ys,fs,step,yoke);
    step = step+1;
end
y = filter(B,1,y);
if step <= 4
    %% MFCC
    overlap = 0.7;

    % Window data
    Window_length_sec = 30/1000; % in seconds
    window_length = fix(Window_length_sec*fs); % calcualting window length
    window_type = hanning(window_length); %type of window function

    s_len = length(ys);

    % In some applications that process large amounts of data with fft, it is common to resize the input so that the number of samples is a power of 2.
    FFT_val = pow2(nextpow2(window_length));
    if FFT_val == window_length
        FFT_length  = pow2(nextpow2(window_length)+1);
    else
        FFT_length = FFT_val;
    end
    %% CEPSTRUM
    %FFT_LEN_V = ceil(FFT_length/2)+1;
    %% Compute Mel filterbank 
    f_lb = 300; f_ub = 4000;
    Mel_low = 1127*log(1+f_lb/700);
    Mel_up = 1127*log(1+f_ub/700);
    Mel_bank = linspace(Mel_low,Mel_up,Mel_coefs+2);
    Hz_bank = 700.*(exp(Mel_bank./1127)-1);
    hhb = Hz_bank;
    Hz_bank = floor((FFT_length +1).*Hz_bank./fs);
    %% Compute Mel filter
    H = mel_bank(FFT_length, Hz_bank);
    cnt = 1;
    %% MFCC gender
    for p = 1:fix(window_length.*(1-overlap)):length(ys)-FFT_length
        signal_fragment = ys(p:p+window_length-1); % getiing part of signal for windowing
        sf_windowed = signal_fragment .* window_type; % applying winow function
        sfw_fft = fft(sf_windowed, FFT_length); % calculate fft of windowed part
        sq_phase = angle(sfw_fft); % Calculate phase of signal
        P_SPECT = abs(sfw_fft).^2/FFT_length;
        %P_SPECT = P_SPECT_A(1:FFT_LEN_V);
        MEL_APP = log(sum(H.*P_SPECT',2)); %% log filterbank energies
        CEPST_COEFS(:,cnt) = dct(MEL_APP);
        %FORMANTS_C(cnt,:) = estimate_formants(sf_windowed,fs);
        cnt = cnt+1;
    end
    %% MFCC age
    if d_on == 1
        cnt = 1;
        for p = 1:fix(window_length.*(1-overlap)):length(y)-FFT_length
            signal_fragment = y(p:p+window_length-1); % getiing part of signal for windowing
            sf_windowed = signal_fragment .* window_type; % applying winow function
            sfw_fft = fft(sf_windowed, FFT_length); % calculate fft of windowed part
            P_SPECT = abs(sfw_fft).^2/FFT_length;
            MEL_APP = log(sum(H.*P_SPECT',2)); %% log filterbank energies
            CEPST_COEFS_2(:,cnt) = dct(MEL_APP);
            cnt = cnt+1;
        end
    end
end
%% Formants
 % Window data
    Window_length_sec = 20/1000; % in seconds
    window_length = fix(Window_length_sec*fs); % calcualting window length
    window_type = hanning(window_length); %type of window function
    preemph      = [1 0.67];
    cnt = 1;
 for p = 1:window_length:length(ys)-window_length
        signal_fragment = ys(p:p+window_length-1); % getiing part of signal for windowing
        sf_windowed = signal_fragment .* window_type; % applying winow function
        fragment = filter(1,preemph,sf_windowed);
        FORMANTS_C(cnt,:) = estimate_formants(fragment,fs,mean(f0));
        cnt = cnt+1;
 end
 
 %%
d_coef = [];
if d_on == 1
    Coefs_delta = CEPST_COEFS_2';%[CEPST_COEFS(1,:); CEPST_COEFS; CEPST_COEFS(end,:)];
    N = 2;
    if d_on
        delta = zeros(size(Coefs_delta,1)-N,Mel_coefs);
        %delta_delta = zeros(size(CEPST_COEFS,1)-2*N,Mel_coefs);
        for d_i = 3:size(Coefs_delta,1)-N
            for n = 1:N
                delta(d_i-2,:) =delta(d_i-2,:)+ n*(Coefs_delta(d_i-n,:)-Coefs_delta(d_i+n,:));
            end 
            delta(d_i-2,:) = delta(d_i-2,:)./(2*(1+4));
        end
%         if size(delta,1) < 2
%             mean_delta = delta;
%         else 
%             mean_delta = mean(delta,1);
%         end
        d_coef = delta;
    else    
    end
end
f0 = mean(f0);
c_coef = CEPST_COEFS';
formants = FORMANTS_C;
end
% if dd_on == 1
%             Coefs_deltadelta = [delta(1,:); delta; delta(end,:)];
%             for dd_i = 3:size(delta)-N
%                 for n = 1:N
%                     delta_delta(dd_i-2,:) =delta_delta(dd_i-n,:)+ n*(Coefs_deltadelta(dd_i-n,:)-Coefs_deltadelta(dd_i+n,:));
%                 end 
%                 delta_delta(dd_i-2,:) = delta_delta(dd_i-2,:)./(2*(1+4));        
%             end
%             if size(delta_delta,1) < 2
%                 mean_deltadelta = delta_delta;
%             else 
%                 mean_deltadelta = mean(delta_delta,1);
%             end
%         end
