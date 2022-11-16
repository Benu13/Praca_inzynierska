function [pi,muu,sigma] = m_step(x,gamma)
%M_STEP Summary of this function goes here
%   Detailed explanation goes here
%%
    N = size(x,1);
    C = size(gamma,2);
    d = size(x,2);
    pi = mean(gamma,1); 
    for i = 1:C
        muu(i,:) = sum(gamma(:,i) .* x,1)/sum(gamma(:,i),1);        
    end
    
    for i = 1:C
        xx = x-muu(i,:);
        sig = sum(gamma(:,i),1)*cov(xx);
        err_flag = 0;
        while err_flag == 0
            try 
                err_flag = 1;
                mvnpdf(x,muu(i,:),sig);
            catch me
                err_flag = 0;
                sig = sig + 0.0001 * eye(rank(x));
                sig = (sig + sig.')/2;
            end
        end

        sigma{i} = sig;
    end
    

end

