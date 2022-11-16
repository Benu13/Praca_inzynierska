function formants = estimate_formants(S3, Fs, t_p)
%%
    % Estymacja formant�w LPC
    % Utworzenie: 18.10.2017, R2017b, J.Przyby�o, AGH
    %


    %  Snell, Roy C., and Fausto Milinazzo. "Formant location from LPC analysis data." 
    %  IEEE� Transactions on Speech and Audio Processing. Vol. 1, Number 2, 1993, pp. 129-134.

    ncoeff=fix(Fs/1000)-10;           % ilo�� wsp�czynnik�w LPC (praktyczna zasada)
    A=lpc(S3,ncoeff);

    % wyznaczenie cz�stotliow�ci formant�w
    r1=roots(A);                  % pierwiastki wielomianu - wyznaczanie
    r=r1(imag(r1)>=0);             % wyb�r rozwi�za� tylko dla zakresu < 0Hz do fs/2
    angz = atan2(imag(r),real(r));
    [frqs,indices] = sort(angz.*(Fs/(2*pi)));
    bw = -1/2*(Fs/(2*pi))*log(abs(r(indices)));

    % wyb�r tylko formant�w o f wi�kszych ni� 90Hz i mniejszych ni� 900Hz
    % https://en.wikipedia.org/wiki/Formant
%     figure()
%     sgtitle('Znajdowanie formant�w z wykorzystaniem LPC');
%     subplot(3,1,1); hold on; grid on;
%     plot(t,S3); title('Sygna� wej�ciowy');
%     xlabel('Czas [s]'); ylabel('Amplituda');
%     subplot(3,1,2); hold on; grid on;
%     fplot(t,A);
%     title('Wielomian predykcji');
%     subplot(3,1,3); hold on; grid on;
%     plot(r1,'o'); plot(r,'*');
    
    formants=[];
    nn = 1;
    for kk = 1:length(frqs)
        if (frqs(kk) > t_p && bw(kk) < 400)
            if kk > 1
                if (frqs(kk)-150) > frqs(kk-1)
                    formants(nn) = frqs(kk);
                    nn = nn+1;
                end
            else
                formants(nn) = frqs(kk);
                nn = nn+1;
            end
        end
    end
    formants = formants(1:3);
end
