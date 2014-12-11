function [res,C,odomlog] = residualRecbot(p,cinfo,mdl,meas,inpt,I,do_dyn)
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

rear_wheel_fi = namesToInds(mdl,{'BL','BR'}); %rear wheel frame indices
steer_fi = namesToInds(mdl,{'steerL','steerR'}); %steer frame indices

% rear_wheel_fi = [6 7];
% steer_fi = [2 3];


ispose = ~isjd; %is body frame pose

if nargin < 7
    do_dyn = false;
end


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
    elseif strcmp(str,'msteer'), m_steer = p(i);
    elseif strcmp(str,'bsteer'), b_steer = p(i);
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

%linear transform for steer
ticksToRad = 3.6864e-6;
steer_mid = m_steer*(inpt.steer*ticksToRad) + b_steer; %steer angle (middle)
speed = inpt.speed;

%calc left/right steer angles for front wheels
fi = steer_fi(1);
L = mdl.frames(fi).HT_parent_jd0(1,4); %assumes no suspension
track_f = 2*mdl.frames(fi).HT_parent_jd0(2,4); %front track width

R = L./tan(steer_mid); %turn radius

steer(:,1) = atan(L./(R-track_f/2)); %left
steer(:,2) = atan(L./(R+track_f/2)); %right
    
%calc left/right velocities for rear wheels (rad/s)
fi = rear_wheel_fi(1);
track_r = 2*mdl.frames(fi).HT_parent_jd0(2,4); %rear track width
rad = mdl.frames(fi).rad; %rear wheel radius

omega = speed./R; %yaw rate

enc(:,1) = (speed - track_r/2*omega)/rad; %left
enc(:,2) = (speed + track_r/2*omega)/rad; %right


%initial & final state
state0 = zeros(ns,1);
state0(isorient) = meas.euler(ind0,:);
state0(ispos) = meas.pos(ind0,:);
state0(steer_fi+dlen) = steer(1,:);

statef = zeros(ns,1);
statef(isorient) = meas.euler(indf,:);
statef(ispos) = meas.pos(indf,:);
statef(steer_fi+dlen) = steer(end,:);

%account for position of pose sensor
state0(ispos) = state0(ispos) - eulerToRot(state0(isorient))*tsensor;

%initial pose covariance
GroundTruthPoseCov = zeros(6); %ground truth pose covariance
P0=zeros(ns);
P0(ispose,ispose) = GroundTruthPoseCov; %covariance of ground truth at time t0


%odometry inputs
odomin = OdometryInput();
odomin.t = time;

%FOR SYSTEMATIC

%controller inputs
odomin.u = NaN(n,nv); %actuated elements of qvel
odomin.u(:,steer_fi+5) = 0;
odomin.u(:,rear_wheel_fi+5) = enc;

vis_act = false(1,nv);
vis_act([steer_fi rear_wheel_fi]+5) = true;
odomin.vis_act = vis_act;
odomin.u = odomin.u(:,vis_act);

%measurements, attitude and steer angles
ismeas = false(1,ns);
odomin.z = NaN(n,ns); %measured elements of state

isatt = false(1,ns); %is attitude
isatt(isorient) = [1 1 0];
odomin.z(:,isatt) = [inpt.rol inpt.pit];
ismeas(isatt) = true;

odomin.z(:,steer_fi+dlen) = steer;
ismeas(steer_fi+dlen) = true;

%FOR STOCHASTIC

%process noise
Qstd = zeros(nv,1); %standard dev. of each element of qvel
Qstd(steer_fi+5) = 1*pi/180;
odomin.Q = diag(Qstd.^2);

%measurement Jacobian
odomin.H = eye(ns);
odomin.H = odomin.H(ismeas,:);

%measurement Covariance
Rstd = NaN(1,ns); %standard deviation
Rstd(isatt) = [1 1]*pi/180;
Rstd(steer_fi+dlen) = .1*pi/180;

odomin.R = diag(Rstd(ismeas).^2);



qvel0 = [];
if do_dyn
    qvel0 = NaN;
end

%do odometry, compute prediction
[statef_pred,odomlog] = odometry(mdl,state0,qvel0,odomin,P0);

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
    h_leg(2) = h;
    set(h_leg(2),'DisplayName','pred')
    set([h; h_end],'Color','g')
    set([he_line; he_cent],'Color','g')
    
    h_leg(1) = plot3(meas.pos(ind0:indf,1), meas.pos(ind0:indf,2), meas.pos(ind0:indf,3), 'b');
    set(h_leg(1),'DisplayName','meas')
    
    legend(h_leg)
    makeLegible(14);
    
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



