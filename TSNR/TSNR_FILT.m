function [sg_filtr_TSNR,sg_filtr_HRNR] = TSNR_FILT(y_s,fs,get_inp,start_end, plot_spec)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%plot_spec = 0;
%addpath('C:/Users/Bennu/Desktop/PR_I_LB/DATA/ROD/SOUND/MONO');
overlap = 0.7;
%signal_file = 'L_a_1.wav' ;
% Window data
Window_length_sec = 30/1000; % in seconds
%FFT_length_times = 2; % length of FFT based on window length (eg. 2 = 2*window_length);

%[y_s, fs] = audioread(signal_file);

%s_len = length(y_s);
window_length = fix(Window_length_sec*fs); % 10ms
window_type = hanning(window_length); %type of window function

% In some applications that process large amounts of data with fft, it is common to resize the input so that the number of samples is a power of 2.
%FFT_length = 2*window_length;
FFT_length  = pow2(nextpow2(window_length));
%% Calculating power spectrum of noise signal
%for every sample in noise signal, we are windowing a part of signal
%starting in that smple and having length of window size. 
%https://www.mathworks.com/help/matlab/math/basic-spectral-analysis.html

start_default = 1;
t_ind = 20*window_length;

if get_inp == 1
    % Getting noise from sample
    dt = 1/fs; % Okres próbkowania 
    ts = 0:dt:(length(y_s)*dt)-dt; % wektor czasu szumu
    figure(1); grid on; hold on;
    if start_end == 1
        title('Select voice starting boundary');
        plot(ts,y_s);
        [ind_x,~]=ginput(2);
        xlabel('Time[s]'); ylabel('Amplitude');
        [~,t_ind] = min(abs(ts-ind_x(2)));
        [~,start_default] = min(abs(ts-ind_x(1)));
        close(figure(1))
    else
        title('Select voice starting boundary');
        plot(ts,y_s);
        [ind_x,~]=ginput(1);
        xlabel('Time[s]'); ylabel('Amplitude');
        [~,t_ind] = min(abs(ts-ind_x));
        close(figure(1))
    end
end
%%

soasv = zeros(FFT_length,1);

for i = start_default:1:t_ind - window_length %length(y_s)-(window_length)
    noise_fragment = y_s(i:i+window_length-1); %selecting noise fragment
    % applaying window to noise fragment
    nf_windowed = noise_fragment .* window_type;
    % calculating sum of absolut squared value of fft 
    soasv = soasv + abs(fft(nf_windowed,FFT_length)).^2;
end
pow_spec = soasv/i;

%PLOT POWER SPECTRUM IF plot_spec == 1
if plot_spec == 1
    figure()
    f = (0:FFT_length-1)*(fs/FFT_length)/10;
    plot(f(1:floor(FFT_length/2)),pow_spec(1:floor(FFT_length/2)))
    xlabel('Frequency')
    ylabel('Power')
end
%% TSNR ALGORITHM 
% POSTERIORI SNR is the ratio of the squared magnitude of the observed noisy signal
% and the noise power. Priori SNR is the ratio of the power of the clean
% signal and of the noise power. For each window in noisy signal we'll
% calculate posteriori SNR and Priori SNR to calculate gain nesessery to
% extract filtered signal from oryginal probe. 
alpha = 0.98; %0.99;
beta = 0.98;
TSNR_OUT = zeros(length(y_s),1);
amp_estimator = zeros(FFT_length,1);
h_val = fix((length(y_s)-FFT_length)/fix(window_length.*(1-overlap)));
sg_amplitude_vec = zeros(FFT_length,h_val);
sq_phase_vec = zeros(FFT_length,h_val);
TSNR_GAIN_VEC = zeros(FFT_length,h_val);
amp_estimator_vec = zeros(FFT_length,h_val);
c = 1;
for p = 1:fix(window_length.*(1-overlap)):length(y_s)-FFT_length % HERE BE DRAGONS
    % DD approach for Wiener filter- needed for TSNR calculations
    signal_fragment = y_s(p:p+window_length-1); % getiing part of signal for windowing
    sf_windowed = signal_fragment .* window_type; % applying winow function
    sfw_fft = fft(sf_windowed, FFT_length); % calculate fft of windowed part
    sg_amplitude = abs(sfw_fft); % Calculate magnitude of signal
    sg_amplitude_vec(:,c) = sg_amplitude;
    sq_phase = angle(sfw_fft); % Calculate phase of signal
    sq_phase_vec(:,c) = sq_phase;
       
    % CALCULATE DD SNR
    post_DD_SNR = sg_amplitude.^2 ./ pow_spec; % calculate posteriori SNR
    prio_DD_SNR = alpha * (abs(amp_estimator.^2) ./ pow_spec) + (1-alpha) * (post_DD_SNR - 1);
    DD_GAIN = prio_DD_SNR ./ (1 + prio_DD_SNR);
    DD_GAIN = max(DD_GAIN,0.1);
    % CALCULATE TSNR 
    prio_TSNR_SNR = abs(DD_GAIN.*sg_amplitude).^2 ./ pow_spec;
    TSNR_GAIN = prio_TSNR_SNR ./ (1+prio_TSNR_SNR);
    TSNR_GAIN_VEC(:,c) = TSNR_GAIN;
    speech_spec_est = TSNR_GAIN .* sg_amplitude;
    amp_estimator = speech_spec_est;
    amp_estimator_vec(:,c) = amp_estimator;
    sfw_fft = amp_estimator.*exp(1i*sq_phase);
    TSNR_OUT(p:p+FFT_length-1) = TSNR_OUT(p:p+FFT_length-1) + real(ifft(sfw_fft,FFT_length))/(1/overlap);
    c = c+1;
end
sg_filtr_TSNR=TSNR_OUT;

S_harmo = max(sg_filtr_TSNR ,0);  %non-linear  function NL(e.g.absolutevalue, minimum or maximum relative to a threshold,etc.)

HRNR_OUT = zeros(length(y_s),1);
c = 1;
for h_loop = 1:fix(window_length.*(1-overlap)):length(y_s)-FFT_length
    % HARMONIC REGENERATION 
    HR_windowed = window_type.*S_harmo(h_loop:h_loop+window_length-1);
    S_harmo_fft = (fft(HR_windowed,FFT_length));
    prio_HRNR_SNR = (TSNR_GAIN_VEC(:,c) .* abs(amp_estimator_vec(:,c)).^2 ...
        + (1 - TSNR_GAIN_VEC(:,c)) .* abs(S_harmo_fft).^2) ./ pow_spec;
    
    HRNR_GAIN = prio_HRNR_SNR ./ (1 + prio_HRNR_SNR);
    HAR_EST = HRNR_GAIN .* sg_amplitude_vec(:,c);
    sfw_fft = HAR_EST.*exp(1i*sq_phase_vec(:,c));
    HRNR_OUT(h_loop:h_loop+FFT_length-1) = HRNR_OUT(h_loop:h_loop+FFT_length-1) + ...
        real(ifft(sfw_fft,FFT_length))/(1/overlap);
    c = c+1;
end
sg_filtr_HRNR=HRNR_OUT;

sg_filtr_TSNR = sg_filtr_TSNR * max(abs(y_s))/max(abs(sg_filtr_TSNR));
sg_filtr_HRNR = sg_filtr_HRNR * max(abs(y_s))/max(abs(sg_filtr_HRNR));

if plot_spec
    figure();
    tiledlayout(1,3)
    nexttile
    [B,f,T] = specgram(y_s,FFT_length,fs,hanning(window_length),window_length-10);
    imagesc(T,f,20*log(abs(B)));axis xy;colorbar
    title(['Spektrogram - mowa zak³ócona'])
    xlabel('Czas [s]');ylabel('Czêstotliwoœæ [Hz]');
    ylim([0 7000])
    xlim([1.5 2]);
    
    %figure();
    nexttile
    [B,f,T] = specgram(sg_filtr_TSNR,FFT_length,fs,hanning(window_length),window_length-10);
    imagesc(T,f,20*log(abs(B)));axis xy;colorbar
    title(['Spektrogram - mowa po filtracji TSNR'])
    xlabel('Czas [s]');ylabel('Czêstotliwoœæ [Hz]');
    ylim([0 7000])
    xlim([1.5 2]);
    
    %figure();
    nexttile
    [B,f,T] = specgram(sg_filtr_HRNR,FFT_length,fs,hanning(window_length),window_length-10);
    imagesc(T,f,20*log(abs(B)));axis xy;colorbar
    title(['Spektrogram - mowa po rekonstrukcji HRNR'])
    xlabel('Czas [s]');ylabel('Czêstotliwoœæ [Hz]');
    ylim([0 7000])
    xlim([1.5 2]);

end
end

