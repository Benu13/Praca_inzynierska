%% Calculating ACF of signal
% Calculate pitch using autocorrelation function


% ACF for long signals is bad so we need to frame signal and then
% calculate ACF for every window 

close all; clear all;
AUDIOPATH = 'C:\Users\Bennu\Desktop\Nowy folder (2)\Semestr_8\Semestr_7\PR_I_LB\DATA\ZDANIA/K1_W_2.wav';
%AUDIOPATH = 'samplen.wav';
[y,fs]=audioread(AUDIOPATH); %Reads .wav file into 2 matrices for amplitude (y) and sampling frequency 
% sprawdzamy czy sygnał jest mono, jeśli nie to zmieniamy
%aa = 'C:\Users\Bennu\Desktop\PR_I_LB\DATA\PL\cv-corpus-5.1-2020-06-22\pl\clips\common_voice_pl_21385179.mp3';
%aa = 'C:\Users\Bennu\Desktop\PR_I_LB\DATA\PL\cv-corpus-5.1-2020-06-22\pl\clips\common_voice_pl_21087650.mp3';
%[y,fs] = audioread(aa);
 
if size(y, 2)>1
    y = mean(y,2); 
end
%detectSpeech(y,fs)
%aa = detectSpeech(y,fs);
%y = y(aa(1,1):aa(1,2));
%% PARAMETERS
freq_range = [80 300]; %frequency of pitch, can help in estimating 
[y,a] = TSNR_FILT(y,fs,1,1);
overlap = 0.5;
% Window data
Window_length_sec = 25/1000; % in seconds

s_len = length(y);
window_length = fix(Window_length_sec*fs); % calcualting window length
window_type = hanning(window_length); %type of window function
N = window_length;
%% MAIN
PC_OUT = zeros(length(y),1);
% P_CEPST_L = zeros(FFT_length,1);
cnt = 1;
PHC_OUT = [];

for p = 1:fix(window_length):length(y)-window_length
    signal_fragment = y(p:p+window_length-1); % getiing part of signal for windowing
    signal_fragment = clipping(signal_fragment);
    sf_windowed = signal_fragment .* window_type; % applying winow function
    preemph = [1 0.64]; sf_filt = filter(1,preemph,sf_windowed); % applying pre emphasis
    R = zeros(N,1);
    for k = 0:N-2
        s_sum = 0;
        for m = 1:N-k-1
            s_sum = s_sum+sf_filt(m)*sf_filt(m+k);
        end
        R(k+1) = s_sum;
    end
    [xx,yy] = findpeaks(R(fix(fs/freq_range(2)):fix(fs/freq_range(1))),'SortStr','descend');
    %% VOICED?UNVOICED
    if length(yy) == 0
        ff0(cnt) = 0;
        VUD(cnt) = 0;
    else
        %ff0(cnt) = fs/((yy(1)+fs/freq_range(2)-1));
        if xx(1) >= 0.55*R(1) & R(1)>=0.1 
            VUD(cnt) = 1;
            ff0(cnt) = fs/((yy(1)+fs/freq_range(2)-1));
            E(cnt) = R(1);
        else
            ff0(cnt) = 0;
            VUD(cnt) = 0;
        end
    end
    cnt = cnt+1;
end
s = ff0(ff0~=0);
mean(s)
%plot(medfilt1(s,5))
%% plots
figure(1)
subplot(2,1,1); hold on; grid on;
title('Sygnał oraz obszar wykrytej mowy');
tt = 0:Window_length_sec:(cnt-2)*Window_length_sec;
tt2 = 0:1/fs:(length(y)/fs)-1/fs;
ix=(VUD~=0);  % this time keep the desired indices
%VUD(VUD == 0) = NaN;
fo = find(ix==1);
% plot speech signal 
plot(tt2,y,'b');
plot(tt,VUD,'r','LineWidth',2) % plot only those particular ones
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
xlim([0 5]);ylim([0 300]);
%%
function s_clipped = clipping(s)
cl = max(s)*0.7;
s_clipped = zeros(length(s),1);
for i = 1:length(s)
    if (s(i) >= cl || s(i) <= -cl)
        s_clipped(i) = s(i);
    else
        s_clipped(i) = 0;
    end
end
end
