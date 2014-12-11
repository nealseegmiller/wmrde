function [res,C,odomlog] = residualSkidsteer(p,cinfo,mdl,meas,inpt,I,do_dyn)
%INPUT
%p:         parameters
%cinfo:     CalibrationInfo object
%mdl:       WmrModel object
%meas:      struct, has .t and other fields
%inpt:      struct, has .t and other fields
%I:         size 2 vector, index in meas.t for t0 and t of interval
%do_dyn:    if true do dynamic sim, else do kinematic sim
%OUTPUT
%res:       nr x 1, residual
%C:         nr x nr, covariance of residual
%odomlog:   OdometryLog object


nf = mdl.nf;
olen=3; %orientation vector length, assumes Euler angles!
dlen=3+olen-1; %ns-nf
ns = nf+dlen; %length of state vector
nv = nf+5;
[isorient,ispos,isjd] = stateIs(ns);

ispose = ~isjd; %is body frame pose



ind0 = I(1); %initial time t0
indf = I(2); %final time t

%interpolate inputs
dt = .04;
[time,inpt] = intervalInputs(meas.t(ind0),meas.t(indf),dt,inpt);
n = length(time);


%PARSE PARAMETERS
cm = mdl.frames(1).cm;

tsensor=zeros(3,1); %position of pose sensor in body coords
%TODO, account for rotation offset of sensor?
yaw0=0; %initial yaw offset


for i = 1:length(cinfo.param_names)
    str=cinfo.param_names{i};
    if strncmp(str,'wgc',3)
        j=eval(str(4:end));
        mdl.wgc_p(j)=p(i);
    elseif strncmp(str,'bsm',3)
        j=eval(str(4:end));
        mdl.bsm_p(j)=p(i);
    elseif strcmp(str,'cmx'), cm(1)=p(i);
    elseif strcmp(str,'cmy'), cm(2)=p(i);
    elseif strcmp(str,'xs'), tsensor(1)=p(i);
    elseif strcmp(str,'ys'), tsensor(2)=p(i);
    elseif strcmp(str,'dyaw0'), yaw0=p(i);
    elseif strncmp(str,'cov',3)
        j=eval(str(4:end));
        mdl.cov_p(j)=p(i);
    end
end

setFrameMass(mdl,1,[],cm,[]);
update_Is(mdl);

meas.euler(:,3) = meas.euler(:,3) + yaw0;


%initial & final state
state0 = zeros(ns,1);
state0(isorient) = meas.euler(ind0,:);
state0(ispos) = meas.pos(ind0,:);

statef = zeros(ns,1);
statef(isorient) = meas.euler(indf,:);
statef(ispos) = meas.pos(indf,:);

%account for position of pose sensor
state0(ispos) = state0(ispos) - eulerToRot(state0(isorient))*tsensor;

%initial pose covariance
GroundTruthPoseCov = zeros(6); %ground truth pose covariance
P0=zeros(ns);
P0(ispose,ispose) = GroundTruthPoseCov; %covariance of ground truth at time t0


%odometry inputs:
odomin = OdometryInput();
odomin.t = time;

%FOR SYSTEMATIC
%controller inputs
odomin.u = inpt.enc;

vis_act = false(1,nv);
vis_act(mdl.wheelframeinds+5) = true;
odomin.vis_act = vis_act;

%measurements
odomin.z = NaN(n,ns); %measured elements of state
isatt = false(1,ns); %is attitude
isatt(isorient) = [1 1 0];
odomin.z(:,isatt) = [inpt.rol inpt.pit];

%FOR STOCHASTIC
%process noise
Qstd = zeros(nv,1); %standard dev. of each element of qvel
odomin.Q = diag(Qstd.^2);

%measurement Jacobian
odomin.H = eye(ns);
odomin.H = odomin.H(isatt,:);

%measurement covariance
Rstd = [1 1]*pi/180; %standard dev. of attitude measurements
odomin.R = diag(Rstd.^2);

if 0
    %DEBUGGING, unicycle
    radii = [mdl.frames(mdl.wheelframeinds).rad];
    rad = radii(1);

    intom = 2.54/100; %inch to meter
    track = (72.5-8)*intom; %track width
    
    Vl = inpt.enc(:,1)*rad;
    Vr = inpt.enc(:,2)*rad;

    Vx = (Vl + Vr)/2;
    Vy = zeros(n,1);
    Wz = (Vr - Vl)/track;

    dVx = 0; dVy = 0; dWz = 0;
    
    c = mdl.bsm_p;    
% %     dVx = c(1)*Vx + c(2)*Wz + c(3)*Vx.*Wz;
%     dVx = c(1)*Vx + c(2)*abs(Wz) + c(3)*Vx.*abs(Wz);
%     dVy = c(4)*Vx + c(5)*Wz + c(6)*Vx.*Wz;
    dWz = c(7)*Vx + c(8)*Wz + c(9)*Vx.*Wz;

    odomin.u(:,4) = Vx + dVx;
    odomin.u(:,5) = Vy + dVy;
    odomin.u(:,3) = Wz + dWz;
    
    %TODO
    [statef_pred,odomlog] = odometryUnicycle(mdl,state0,odomin,P0);
    
else
    
    qvel0 = [];
    if do_dyn
        qvel0 = NaN;
    end

    [statef_pred,odomlog] = odometry(mdl,state0,qvel0,odomin,P0);
end


%account for position of pose sensor
statef_pred(ispos) = statef_pred(ispos) + eulerToRot(statef_pred(isorient))*tsensor;

set_HT_sensor_to_body(odomlog,poseToHT(zeros(3,1),tsensor));

if 0
    %%
    %DEBUGGING
    I = floor(linspace(1,n,5));
    do_sensor = true;
    dim = 3;
    marg = [];
    [h,h_end,he_line,~,he_cent] = plotPath(odomlog,do_sensor,I,dim,marg);
    set([h; h_end; he_line; he_cent],'Color','g')
    set(he_line,'LineStyle',':')
    plot3(meas.pos(ind0:indf,1),meas.pos(ind0:indf,2),meas.pos(ind0:indf,3),'b')
    
    
    set(figure,'name','attitude')
    hold on
    plot(meas.t(ind0:indf),meas.euler(ind0:indf,[1 2])*180/pi);
    plot(time,odomlog.state(:,[1 2])*180/pi,':')
    xlabel('time (s)')
    ylabel('deg')
    legend('meas rol','meas pit','est rol','est pit','Location','Best')
    makeLegible(14)
    
end



%compute residual: meas - pred
inres = cinfo.inresidual; %abbreviate
res = statef(inres) - statef_pred(inres);  %requires angles not wrapped!

%compute residual covariance
C = odomlog.P(:,:,end);
C(ispose,ispose) = C(ispose,ispose) + GroundTruthPoseCov; %covariance of ground truth at time t
C = C(inres,inres);

if any(cinfo.isstoch)
    %do stochastic calibration, append stochastic residual & covariance
    [res,C]=appendStoch(res,C);
end



