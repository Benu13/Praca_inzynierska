[y_s,fs] = audioread('samplen.wav');
[a,b] = TSNR_FILT(y_s(:,1),fs,1,0,1);
T = 1/fs;
t = 0:T:(length(y_s)*T)-T;

figure('Name','Przebieg sygna�u w czasie.'); hold on;
sgtitle('           Przebieg sygna��w przed i po filtracji:');

subplot(3,1,1); hold on; grid on;
plot(t, y_s); 
title('Sygna� zak��cony:');
xlabel('Czas [s]'); ylabel('Amplituda');

subplot(3,1,2); hold on; grid on;
plot(t, a); 
title('Filtracja przy wykorzystaniu metody TSNR:');
xlabel('Czas [s]'); ylabel('Amplituda');

subplot(3,1,3); hold on; grid on;
plot(t, b); 
title('Sygna� po rekonstrukcji HRNR:');
xlabel('Czas [s]'); ylabel('Amplituda');