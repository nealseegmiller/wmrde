function [x,P,md] = kfMeasUpdate(x,P,y,H,R,md_tol)
%INPUTS
%x:         (n x 1) state estimate
%P:         (n x n) state cov
%y:         (m x 1) innovation/residual
%H:         (m x n) Jacobian
%R:         (m x m) meas cov
%md_tol:    Mahalanobis Distance threshold, no update if MD > thresh
%           if < 1 md_tol is a probability level, thresh = chi2inv(1-plevel,length(y))
%           e.g. for plevel = .001, and 2 degrees of freedom, thresh = 13.82


if nargin < 6
    md_tol = Inf;
end

if md_tol < 1
    thresh = chi2inv(md_tol,length(y));
else
    thresh = md_tol;
end
    
S = H*P*H' + R; %residual covariance
K = P*H'/S; %Kalman gain
md = y'*(S\y); %Mahalanobis distance
if md < thresh;
    x = x + K*y;
    P = (eye(length(x))-K*H)*P;
end