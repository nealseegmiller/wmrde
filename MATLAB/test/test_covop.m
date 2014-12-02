% test_covop
%test function for computing covariance of outer product of random vector
clear_all
close all

setseed(852)


%symmetric random matrix
d=3; %length of random vector
c = rand(d);
c = c*c';

C=covop(c)

%verify numerically

n = 1e6; %number of samples
%sample
s = mvnrnd(zeros(n,d),c);
%compute outer product
op = zeros(n,(d^2-d)/2+d);
for i = 1:n
    op(i,:) = triuv(s(i,:)'*s(i,:));
end

cov(op)
C-cov(op)

