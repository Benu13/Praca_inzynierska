function [ff0, v_comb] = ACF_FF(y,fs)
%% Calculating ACF of signal 
% ACF for long signals is bad m'key so we need to frame signal and then
% calculate ACF for every window 
%% PARAMETERS
freq_range = [100 280]; %frequency of pitch, can help in estimating 

% Window data
Window_length_sec = 40/1000; % in miliseconds

s_len = length(y);
window_length = fix(Window_length_sec*fs); % calcualting window length
window_type = hanning(window_length); %type of window function
N = window_length;
%% MAIN
cnt = 1;
ff0 = zeros(fix(s_len/window_length),1);
v_comb = [];
for p = 1:fix(window_length):s_len-window_length
    signal_fragment = y(p:p+window_length-1); % getiing part of signal for windowing
    signal_fragment_c = clipping(signal_fragment);
    sf_windowed = signal_fragment_c .* window_type; % applying winow function
    sf_filt = sf_windowed;
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
    if isempty(yy)
        ff0(cnt) = 0;
    else
        if xx(1) >= 0.45*R(1) && R(1)>=0.1
            ff0(cnt) = fs/((yy(1)+fs/fix(freq_range(2))-1));
            v_comb = [v_comb; signal_fragment];
        end
    end
    cnt = cnt+1;
end
ff0 = ff0(ff0~=0);
%%
function ss = clipping(s)
cl = max(s)*0.6;
ss = zeros(length(s),1);
for i = 1:length(s)
    if (s(i) >= cl || s(i) <= -cl)
        ss(i) = s(i) - sign(ss(i))*cl ;
    else
        ss(i) = 0;
    end
end
end

end

