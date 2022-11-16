%% Calculating CEPSTRUM of signal  
 clear all;
AUDIOPATH = 'C:/Users/Bennu/Desktop/PR_I_LB/DATA/Zdania/M1_W_22.wav';
%AUDIOPATH = 'C:/Users/Bennu/Desktop/PR_I_LB/DATA/SGL/K1_a_2.wav';
%AUDIOPATH = 'samplen.wav';

[y,fs]=audioread(AUDIOPATH); %Reads .wav file into 2 matrices for amplitude (y) and sampling frequency 
% sprawdzamy czy sygnał jest mono, jeśli nie to zmieniamy
if size(y, 2)>1
    y = mean(y,2); 
end
[~,y] = TSNR_FILT(y,fs,1,0);
%% PARAMETERS
freq_range = [80 160]; %frequency of pitch, needed for component at zero quefrency filtration 
quefency = fix([fs/freq_range(2) fs/freq_range(1)]);
%[y,a] = TSNR_FILT(y,fs,1,1);
overlap = 0;
% Window data
Window_length_sec = 40/1000; % in seconds

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
    E(cnt) = xx./P_CEPST(1);
    if xx/P_CEPST(1)>= 0.00015
        ff0(cnt) = fs/(yy+quefency(1));
        VUD(cnt) = 1;
    else
        ff0(cnt) = 0;
        VUD(cnt) = 0;
    end
    cnt = cnt+1;
end

figure()
subplot(2,1,1); hold on; grid on;
title('Sygnał oraz obszar wykrytej mowy');
tt = 0:Window_length_sec:(cnt-2)*Window_length_sec;
tt2 = 0:1/fs:(length(y)/fs)-1/fs;
ix=(VUD~=0);  % this time keep the desired indices
%VUD(VUD == 0) = NaN;
fo = find(ix==1);
% plot speech signal 
plot(tt2,y,'b');
for i = 2:length(VUD)-1
    VUD_ol = VUD(i-1);
    VUD_fut = VUD(i+1);
    if VUD(i) == 1
        plot([tt(i)-Window_length_sec/2 tt(i)+Window_length_sec/2], [1 1],'r','LineWidth',2);
        plot([tt(i)-Window_length_sec/2 tt(i)+Window_length_sec/2], [-1 -1],'r','LineWidth',2);
    end
    if ~VUD_ol && VUD(i) 
        xline(tt(i)-Window_length_sec/2,'r','LineWidth',2);
    end
    if ~VUD_fut && VUD(i) 
        xline(tt(i)+Window_length_sec/2,'r','LineWidth',2);
    end
    
    if i==2 && VUD_ol
        xline(tt(i-1),'r','LineWidth',2);
    end
end
xlim([0 tt(end)]);
% plot(tt,VUD,'r','LineWidth',2) % plot only those particular ones
legend('Przebieg sygnału','Obszar wykrytej mowy');
xlabel('Czas [s]');ylabel('Amplituda');
%plot([tt(fo(1)) tt(fo(1))],[0 1],'r','LineWidth',2,'HandleVisibility','off')
%plot([tt(fo(end)) tt(fo(end))],[0 1],'r','LineWidth',2,'HandleVisibility','off')
subplot(2,1,2);hold on; grid on;
title('Wykres tonu podstawowego');
xlabel('Czas [s]'); ylabel('Częstotliwość [Hz]');
ff0(ff0==0) = NaN;
plot(tt,ff0,'b')
scatter(tt(ix),ff0(ix),'b*')
xlim([0 tt(end)]);ylim([0 300]);
%%
figure(3)
sgtitle('Stosunek gamnitude dla quefrency tonu podstawowego do gamnitude dla quefrency 0');
subplot(2,2,2); hold on; grid on;
title("Wykres czasowy samogłoski ""a"" dla kobiety");
tt = 0:Window_length_sec:(cnt-2)*Window_length_sec;
plot(tt2,y,'b');
xlabel('Czas [s]');ylabel('Amplituda');

subplot(2,2,4); hold on; grid on;
title('Otrzymany stosunek');
tt = 0:Window_length_sec:(cnt-2)*Window_length_sec;
tt2 = 0:1/fs:(length(y)/fs)-1/fs;
plot(tt,E,'b');
yline(4*10^-4,'r');
xlabel('Czas [s]');ylabel('Wartość');
legend('Stosunek', 'P_v_u = 4*10^-^4');
%%
% aa = 0;
% for i = 1:cnt-1
%     hold on;
%     ss = locks{i,:};
%     scatter(i,fs/ss(4))
%     aa = aa + ss(4);
% end
% %% LIFTERING
%xx = estimate_formants(y,fs);