function loss = cl_fun(x,pi,muu,sigma,gamma)
%CL_FUN Summary of this function goes here
%   Detailed explanation goes here
%%
N = size(x,1);
C= size(gamma,2);
loss =[];
for i = 1:C
    dista = mvnpdf(x,muu(i,:),cell2mat(sigma(i)));
    %dista = multi_dist(x,pi(i),muu(i,:),cell2mat(sigma(i)));
    loss(:,i) = log(pi(i)*dista);
%     dista = log(mvnpdf(x,mu(i,:),cell2mat(sigma(i))));
%     loss(:,i) = gamma(:,i) * log(pi(i)+0.00001)+dista-log(gamma(:,i)+0.000001);
end

loss = sum(loss,'All');
end
%    dista = mvnpdf(x,mu(i,:),cell2mat(sigma(i)));
  %  loss(:,i) = log(pii(i)*dista);