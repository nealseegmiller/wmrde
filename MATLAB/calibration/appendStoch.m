function [r,R] = appendStoch(r,R,instoch)
%append stochastic residual and measurement covariance to systematic
%INPUTS
%r:         (nr x 1) systematic residual
%R:         (nr x nr) systematic measurement covariance
%in:        (optional) size nr logical, which elements of r to include in stochastic residual
%OUTPUTS
%r,R:       with stochastic values appended

if nargin < 3
    instoch = true(size(r));
end

nr = length(r);

r_ = r(instoch);
R_ = R(instoch,instoch);

rstoch = triuv(r_*r_') - triuv(R_); %meas - pred
Rstoch = covop(R_);

%append
r=[r; rstoch];
Rsys=R; %backup
R = zeros(length(r)); %block diagonal
R(1:nr,1:nr) = Rsys;
R(nr+1:end,nr+1:end) = Rstoch;