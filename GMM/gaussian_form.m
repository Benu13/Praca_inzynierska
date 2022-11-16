function gf = gaussian_form(x,p,mean_vec, cov_mat)
%GAUSSIAN_FORM Summary of this function goes here
%   Detailed explanation goes here
%%
mean_vec = mu(i,:);
cov_mat = cell2mat(sigma(i));
p = pi(i);
gf = (1./sqrt(((2*p).^ size(x,2))*det(cov_mat))) .* exp(-0.5.*...
    (x-mean_vec) * inv(cov_mat) * (x-mean_vec)');
end

