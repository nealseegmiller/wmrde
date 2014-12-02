function [state,odomlog] = odometry(mdl,state0,qvel0,odomin,P0)
%odometry
%INPUTS
%state0:    ns x 1, state at initial time t0
%qvel0:     nv x 1, joint space velocity at t0, empty for kinematic sim, NaN for automatic calc
%odomin:    OdometryInput object
%P0:        (optional) ns x ns, state covariance at time t0
%OUTPUTS
%state:     ns x 1, state at final time t
%odomlog:   OdometryLog object

if nargin < 5
    P0 = [];
end

if isempty(qvel0)
    do_dyn = false; %do kinematic sim
else
    do_dyn = true; %do dynamic sim
    
    qvel = qvel0;
    
    mdl.use_constraints = 1;
    mdl.act_fh = []; %ideal actuators
end

if isempty(P0)
    do_stoch = false;
else
    do_stoch = true;
    P = P0;
end


%for indexing
nf = mdl.nf;
nw = mdl.nw;
olen=3; %orientation length
dlen=3+olen-1;
ns = nf+dlen; %number of states
nv = nf-1+6;
isorient = stateIs(ns,nf);

wsi = mdl.wheelframeinds+dlen; %indices of wheel rotations in state vector

nsteps = length(odomin.t);

radii = [mdl.frames(mdl.wheelframeinds).rad];

%initialize
state = state0;
dz = mdl.dz_target;


odomlog=[];
if nargout > 1
    %DATA LOGGING
    odomlog = OdometryLog(nsteps,ns,nv,nw,do_stoch);
    set_time(odomlog,1:nsteps,odomin.t)
    set_state(odomlog,1,state);
    if do_stoch
        set_P(odomlog,1,P)
    end
end


iter = 0; %forwardDyn will overwrite

for i = 1:(nsteps-1)
       
    %compute transforms, update wheel/terrain contact frames
    [HT_parent, HT_world] = stateToHT(mdl,state);

    %compute contact angles
    if strcmp(mdl.description,'crusher')
        contact_angles = crusherContactAngles(mdl,state,HT_world); %DEBUGGING
    else
        contact_angles = -state(wsi)'; %contact_angles_rot0=0
    end
    
    contacts = updateWheelContactGeomBlind(radii, HT_world(:,:,mdl.wheelframeinds), contact_angles);
    %copy dz
    for wno = 1:nw
        contacts(wno).dz = dz(wno);
    end
        
    dt = odomin.t(i+1)-odomin.t(i);
    u = odomin.u(i,:)';

    if do_dyn
        
        if isnan(qvel)
            %init qvel using kinematic model
            qvel = forwardVelKin(mdl,state,u,HT_world,contacts);
        end
        
        [qacc,~,fdout] = forwardDyn(mdl, state, qvel, u, [], HT_parent, HT_world, contacts, dt);
        qvel = qvel + qacc*dt; %symplectic Euler
        
        vc = fdout.vc0; %???
    else
        
        [qvel,vc] = forwardVelKin(mdl,state,u,HT_world,contacts); %compute velocity, DAE
    end
    %integrate velocity
    statedot = qvelToQdot(qvel, state(isorient), HT_world(1:3,1:3,1));
    state = state + statedot*dt;
    
    if do_dyn %TODO, not for kinematic?
        
        dz = dz + vc(3,:)*dt; %update dz
        dz = dz + (dt/.1)*(mdl.dz_target - dz); %DEBUGGING, for drift
    end
    

    z = odomin.z(i+1,:)';
    ismeas = ~isnan(z);
    z = z(ismeas);
    
    if do_stoch
        
        %predict
        Q_ = odomin.Q;
        Q_(1:6,1:6) = Q_(1:6,1:6) + bodyVelNoise(mdl.cov_p,qvel(1:6)); %add noise to body frame velocity
        P = predictStateCov(P,Q_,state(isorient),HT_world(1:3,1:3,1),qvel(1:6),dt);
        
        %measurement update
        y = z - state(ismeas); %innovation
        H_ = odomin.H;
        R_ = odomin.R;
        
        [state,P] = kfMeasUpdate(state,P,y,H_,R_,Inf);
        
    else
        state(ismeas) = z;
        
    end


    if ~isempty(odomlog)
        %DATA LOGGING
        set_state(odomlog,i+1,state);
        set_qvel(odomlog,i+1,qvel);
        %DEBUGGING
        set_vc(odomlog,i,vc);
        set_contacts(odomlog,i,contacts);
        set_iter(odomlog,i,iter);
        if do_stoch
            set_P(odomlog,i+1,P);
            set_Q(odomlog,i+1,Q_);
        end
    end

    
end




