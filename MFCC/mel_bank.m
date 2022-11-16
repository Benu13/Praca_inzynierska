function MelF_co = mel_bank(F_len,mel_bank)
%MEL_BANK - create mel-frequency bank
MelF_co = zeros(length(mel_bank)-2,F_len);
for m = 2:length(mel_bank)-1
    for k = 1:F_len-1
        if k<mel_bank(m-1)
            MelF_co(m-1,k) = 0;
        elseif (k>=mel_bank(m-1)&&k<=mel_bank(m))
            MelF_co(m-1,k) = (k-mel_bank(m-1))/(mel_bank(m)-mel_bank(m-1));
        elseif (k>=mel_bank(m)&&k<=mel_bank(m+1))
            MelF_co(m-1,k) = (mel_bank(m+1)-k)/(mel_bank(m+1)-mel_bank(m));
        elseif k>mel_bank(m+1)
            MelF_co(m-1,k) = 0;    
        end
    end
end


end

