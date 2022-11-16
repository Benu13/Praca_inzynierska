function gamma = e_step(x,c,pi, muu, sigma)
%E_STEP Summary of this function goes here
%   Detailed explanation goes here
N = size(x,1);
for i = 1:c
    gamma(:,i) = pi(i) * mvnpdf(x,muu(i,:),cell2mat(sigma(i)));
    %gamma(:,i) = pi(i) * multi_dist(x,pi(i),muu(i,:),cell2mat(sigma(i)));
end

gamma = gamma ./ sum(gamma,1);
end

