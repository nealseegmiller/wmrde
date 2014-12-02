function [cinfo,x,P,Q] = initcalSkidsteer(mdl,do_dyn)
%INPUT
%mdl:       WmrModel object
%do_dyn:    dynamic simulation (else kinematic)
%OUTPUT
%cinfo:     CalibrationInfo object
%x:         np x 1, initial state in identification filter (i.e. parameter values)
%P:         np x np, initial state covariance
%Q:         np x np, process noise

cinfo = CalibrationInfo();

%SYSTEMATIC
if do_dyn
    nsp = numel(mdl.wgc_p);
    param_names_sys = catstrnum('wgc',nsp);
    x_sys = mdl.wgc_p(:)';
    Pstd_sys = abs(mdl.wgc_p')/2; %standard deviations
else
    %body slip model
    nsp=9; %number of body slip model parameters
    param_names_sys = catstrnum('bsm',nsp);
    x_sys = zeros(1,nsp);
    Pstd_sys = .1*ones(1,nsp);
%     Pstd_sys = [.1 .1 .1 .01 .01 .05]; %DEBUGGING, smaller for angular slip params
end
%append params for sensor location, yaw initialization error
param_names_sys = [param_names_sys,'cmx','cmy','xs','ys','dyaw0'];
x_sys = [x_sys, 0, 0, 0, 0, 0];
Pstd_sys = [Pstd_sys, .1, .1, .2, .2, 2*pi/180];

nsys=length(x_sys);

%STOCHASTIC 
x_stoch = [.2 .2 [4 2]*pi/180]; %LandTamer
% x_stoch = [.1 .3 [1.5 2]*pi/180]; %DEBUGGING, LandTamer on parkinglot
% x_stoch = [.1 .5 [7 2]*pi/180]; %DEBUGGING, LandTamer on dirt?
% x_stoch = [1.2 .7 [10 20]*pi/180]; %DEBUGGING, Crusher at Camp Roberts, md thresh = 90


%IMPORTANT THAT ATTITUDE NOISE IS SUFFICIENTLY LARGE! (RELATIVE TO MEASUREMENT COVARIANCE)
%OTHERWISE ATTITUDE MEASUREMENT UPDATES WILL AFFECT POSITION

nstoch=length(x_stoch);

cinfo.param_names = [param_names_sys, catstrnum('cov',nstoch)];
cinfo.isstoch = [false(1,nsys), true(1,nstoch)]; %indicates if parameter is stochastic
x = [x_sys,x_stoch]';
Pstd = [Pstd_sys, x_stoch/2];
Qstd = [zeros(1,nsys), 1e-6*ones(1,nstoch)]; %per meter

P = diag(Pstd.^2);
Q = diag(Qstd.^2);


%EDIT THIS

cinfo.isfixed = true(1,nsys+nstoch);
if do_dyn
    %TODO
    cinfo.isfixed(5:7) = false; %odeWgc, Cx,Cy,Crr
    cinfo.isfixed(nsp+(1:4)) = false; %cmx cmy xs ys
    
else
%     cinfo.isfixed(1:nsp) = false;
%     cinfo.isfixed(nsp+[3 4]) = false; %xs, ys
%     cinfo.isfixed(8) = false; %DEBUGGING

end
% cinfo.isfixed(cinfo.isstoch) = [false(1,3) true];


%logical, indicates which elements of state are in calibration residual
dlen=6-1;
ns=mdl.nf+dlen; %number of states)
cinfo.inresidual=false(1,ns);
[isorient,ispos] = stateIs(ns,mdl.nf);
cinfo.inresidual(isorient) = [0 0 1]; %yaw
cinfo.inresidual(ispos) = [1 1 0]; %x,y





