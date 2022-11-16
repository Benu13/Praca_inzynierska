function y = multi_dist(x,pi,muu,sigma)
%MULTI_DIST Summary of this function goes here
%   Detailed explanation goes here
for i = 1:size(x,1)
    y(i,:) = (1./sqrt(((2*pi).^size(x,2))*det(sigma))) .* exp(-0.5.*...
    (x(i,:)-muu) * inv(sigma) * (x(i,:)-muu)');
end
end

