%% MFCC, Mel-frequency cepstral coefficients
% Najpierw obliczymy ton podstawowy z wykorzystaniem ACF, przy okazji
% wykorzystamy decyzję czy sygnał jest voiced/unvoiced żeby stworzyć
% kombinację samych ramek które zostały określone jako voiced co pozwoli
% nam pominąć problem szumu. Na poswstałym zlepku będziemy badać MFCC. 
%% Signal reading
clear all;
AUDIOPATH = 'C:/Users/Bennu/Desktop/PR_I_LB/DATA/SGL/K1_a_1.wav';
%AUDIOPATH = 'samplen.wav';
[y,fs]=audioread(AUDIOPATH); %Reads .wav file into 2 matrices for amplitude (y) and sampling frequency 
% sprawdzamy czy sygnał jest mono, jeśli nie to zmieniamy
if size(y, 2)>1
    y = mean(y,2); 
end
%% TSNR_filtration
% [a,y] = TSNR_FILT(y,fs);

%% Get pitch
[f0, y] = ACF_FF(y,fs,1,0);
aff=y;
%%
B = [1 -0.98];
y = filter(B,1,y);
%% MFCC
overlap = 0.5;

% Window data
Window_length_sec = 25/1000; % in seconds
window_length = fix(Window_length_sec*fs); % calcualting window length
window_type = hanning(window_length); %type of window function

s_len = length(y);

% In some applications that process large amounts of data with fft, it is common to resize the input so that the number of samples is a power of 2.
FFT_val = pow2(nextpow2(window_length));
if FFT_val == window_length
    FFT_length  = pow2(nextpow2(window_length)+1);
else
    FFT_length = FFT_val;
end
%% CEPSTRUM

PC_OUT = zeros(length(y),1);
rec = zeros(length(y),1);
FFT_LEN_V = ceil(FFT_length/2)+1;
cnt = 1;
%% Compute Mel filterbank 
Mel_coefs = 10;
f_lb = 300; f_ub = 6000;
Mel_low = 1127*log(1+f_lb/700);
Mel_up = 1127*log(1+f_ub/700);
Mel_bank = linspace(Mel_low,Mel_up,Mel_coefs+2);
Hz_bank = 700.*(exp(Mel_bank./1125)-1);
Hz_bank = floor((FFT_length +1).*Hz_bank./fs);
%% Compute Mel filter
H = mel_bank(FFT_length, Hz_bank);
cnt = 1;
%% MFCC
for p = 1:fix(window_length.*(1-overlap)):length(y)-FFT_length
    signal_fragment = y(p:p+window_length-1); % getiing part of signal for windowing
    %signal_fragment = filter(B,1,signal_fragment);
    signal_fragment2 = aff(p:p+window_length-1); % getiing part of signal for windowing
    sf_windowed = signal_fragment .* window_type; % applying winow function
    sf_windowed2 = signal_fragment2 .* window_type; % applying winow function
    sfw_fft = fft(sf_windowed, FFT_length); % calculate fft of windowed part
    sfw_fft2 = fft(sf_windowed2, FFT_length); % calculate fft of windowed part
    
    sq_phase = angle(sfw_fft); % Calculate phase of signal
    P_SPECT = abs(sfw_fft).^2/FFT_length;
    P_SPECT2 = abs(sfw_fft2).^2/FFT_length;
    %P_SPECT = P_SPECT_A(1:FFT_LEN_V);
    MEL_APP = log(sum(H.*P_SPECT',2)); %% log filterbank energies
    CEPST_COEFS(:,cnt) = dct(MEL_APP);
    cnt = cnt+1;
    if cnt == 2
        ssf = signal_fragment;
        fssw = sf_windowed;
        aaa = sfw_fft;
        ppps = P_SPECT;
        mapp = CEPST_COEFS(:,cnt-1);
        kek = P_SPECT2;
    end
end
%% Plot 3d MFCC work
% Plot signal
figure()
n = FFT_length;
f = (0:n-1)*(fs/n);
t = (0:length(ssf)-1)/fs;
tiledlayout(2,2);

nexttile; hold on; grid on;
plot(t,ssf);
title('Pobrany fragment sygnału');
xlabel('Czas [s]'); ylabel('Amplituda');
xlim([-2/1000 27/1000]); ylim([-0.3 0.3]);

nexttile; hold on; grid on;
plot(t,fssw)
xlabel('Czas [s]'); ylabel('Amplituda');
xlim([-2/1000 27/1000]); ylim([-0.20 0.20]);
title('Okienkowanie fragmentu z wykorzystaniem funkcji Hanna');

nexttile; hold on; grid on;
for i= 1:size(H,1)
    plot(f(1:floor(n/5)),(H(i,1:floor(n/5))))
end
xlabel('Częstotliwość [Hz]'); ylabel('Amplituda');
title('Przebieg filtrów Melowych');

nexttile; hold on; grid on;
plot(f(1:floor(n/5)),mag2db(ppps(1:floor(n/5))));
xlabel('Częstotliwość [Hz]'); ylabel('Gęstość widmowa mocy [dB/Hz]');
title('Gęstość widmowa mocy próbki');
%%
fig = figure;

hold on; grid on;
%set(gca,'xtick',[0:1:Mel_coefs+1]);
xlim([0 11]);
xlabel('Numer filta'); ylabel('Częstotliwość [Hz]'); zlabel('Gęstość widmowa mocy [dB/Hz]');
title('Wynik filtracji PCB filtrem Mel');

po = H.*P_SPECT';
cnt = 1;
for i = 1:2:2*5
    subplot(5,2,i); hold on; grid on;
    %plot3(i*ones(floor(n/5),1),f(1:floor(n/5)),mag2db(po(i,1:floor(n/5))))
    plot(f(1:floor(n/5)),mag2db(po(cnt,1:floor(n/5))))
    title(['Filtr nr. ' num2str(cnt)]);
    cnt = cnt+1;
end
xlabel('Częstotliwość [Hz]','FontSize',16,'FontWeight','normal');
subplot(5,2,[1:5]*2)
hold on; grid on;
bar(1:5,mapp(1:5));
title('Wartości współczynników cepstralnych');
ylabel('Wartość współczynnika','FontSize',12,'FontWeight','normal');
xlabel('Numer współczynnika','FontSize',12,'FontWeight','normal');
set(gca,'ytick',[-14:1:4]);
hold off

han=axes(fig,'visible','off'); 
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Gęstość widmowa mocy [dB/Hz]','FontSize',16,'FontWeight','normal');
sgtitle('Wynik filtracji oraz obliczone wartości cepstralne','FontSize',16,'FontWeight','bold');

%%
figure()
sgtitle('Gęstość widmowa mocy próbek');
subplot(1,2,2);%hold on; grid on;
%plot(f(1:floor(n/5)),mag2db(ppps(1:floor(n/5))));
semilogx(f(1:floor(n/5)),mag2db(ppps(1:floor(n/5))));
xlabel('Częstotliwość [Hz]'); ylabel('Gęstość widmowa mocy [dB/Hz]');
title('Próbka po zastosowaniu filtra preemfazy');grid on;
% ylim([-180 20]);
subplot(1,2,1);%hold on; grid on;
%plot(f(1:floor(n/5)),mag2db(kek(1:floor(n/5))));
semilogx(f(1:floor(n/5)),mag2db(kek(1:floor(n/5))));
xlabel('Częstotliwość [Hz]'); ylabel('Gęstość widmowa mocy [dB/Hz]');
title('Próbka przed zastosowaniem filtra preemfazy');grid on;
% ylim([-180 20]);