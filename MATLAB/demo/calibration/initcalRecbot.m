function [cinfo,x,P,Q] = initcalRecbot(mdl,do_dyn)
%INPUTS
%mdl:       WmrModel object
%do_dyn:    dynamic simulation (else kinematic)

cinfo = CalibrationInfo();

%SYSTEMATIC
% dyn=true; %DEBUGGING
if do_dyn
    %handle multi-column wgc_p
    nsp = numel(mdl.wgc_p);
    param_names_sys = catstrnum('wgc',nsp);
    x_sys = mdl.wgc_p(:)';
    Pstd_sys = abs(x_sys)/2;
else
    %body slip model
    nsp=6; %number of slip model parameters
    param_names_sys = catstrnum('bsm',nsp);
    x_sys = zeros(1,nsp);
    Pstd_sys = .1*ones(1,nsp);
end
%append params for sensor location, yaw initialization error
param_names_sys = [param_names_sys,'msteer','bsteer','cmx','cmy','xs','ys','dyaw0'];
x_sys = [x_sys, 1, 0, 0, 0, 0, 0, 0];
if do_dyn %cmx
    %TODO, for kinematic too?
    x_sys(nsp+3) = 1.770/2;
end

Pstd_sys = [Pstd_sys, .1, .1, .1, .1, .2, .2, 2*180/pi];


nsys=length(x_sys);

%STOCHASTIC 
x_stoch = [.3 .05 [4 10]*pi/180];  %TODO

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
    %TODO, why does ode_wgc do better than pacjeka wgc?
    cinfo.isfixed(5:7) = false; %ode_wgc Cx,Cy,Crr
    cinfo.isfixed(nsp+(1:6)) = false; %msteer bsteer cmx cmy xs ys
    
%     cinfo.isfixed([4 9 14]) = false; %pacejka_wgc Bx,By,Crr
else
%     cinfo.isfixed(1:nsp) = false;
%     cinfo.isfixed(nsp+[1 2 5 6]) = false; %msteer bsteer xs ys
end
% cinfo.isfixed(cinfo.isstoch) = [false(1,3) true];


%logical, indicates which elements of state are in calibration residual
dlen=6-1;
ns=mdl.nf+dlen; %number of states)
cinfo.inresidual=false(1,ns);
[isorient,ispos] = stateIs(ns,mdl.nf);
cinfo.inresidual(isorient) = [0 0 1]; %yaw
cinfo.inresidual(ispos) = [1 1 0]; %x,y





