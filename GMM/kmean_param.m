function [initial_pi,initial_means,initial_cov] = kmean_param(x,clusters)
%KMEAN_PARAM Summary of this function goes here
%   Detailed explanation goes here
predictions = kmeans(x,clusters,'MaxIter',1000);
% init_means = zeros(
%%
labels = unique(predictions);
for l = 1:length(labels)
    idx = find(predictions == labels(l));
    initial_pi(l,:) = size(idx,1)./ size(x,1);
    initial_means(l,:) = mean(x(idx,:),1);
    d_meaned = x(idx,:) - initial_means(l,:);
    Nk = size(x(idx),1);
    sigma = ((initial_pi(l,:) *d_meaned.') * d_meaned)/ Nk;
    
    err_flag = 0;
    while err_flag == 0
        try 
            err_flag = 1;
            mvnpdf(x,initial_means(l,:),sigma);
        catch me
            err_flag = 0;
            sigma = sigma + 0.0001 * eye(rank(x));
            sigma = (sigma+sigma')/2;
        end
    end
        
    initial_cov{l} =sigma;
end

end

