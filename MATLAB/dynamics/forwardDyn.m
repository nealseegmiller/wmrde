function [qacc, u_err, out] = forwardDyn(mdl, state0, qvel0, u_cmd, u_interr, HT_parent, HT_world, contacts, dt)
%forward dynamics, compute joint space acceleration.
%use suffix 0 in variable names to distinguish values at time i from values at time i+1, as needed
%INPUTS
%mdl:       WmrModel object
%state0:    ns x 1, required only as input for joint constraint function (assumed to be consistent with HT_parent, HT_world)
%qvel0:     nv x 1, joint space velocity
%u_cmd:     nv x 1, commanded velocities for actuator model, NaN for not actuated DOF
%u_interr:  na x 1, integrated error for actuator model
%HT_parent: 4x4x nf, Homogeneous transforms to parent coords, HT_parent(:,:,i) = HT^(parent)_(frame i)
%HT_world:  4x4x nf, to world coords, HT_world(:,:,i) = HT^(world)_(frame i)
%contacts:  1 x nw array of WheelContactGeom objects, or
%           1 x nt array of TrackContactGeom objects
%dt:        integration step size (only needed if using constraints)
%OUTPUTS
%qacc:      ns x 1, joint space acceleration, d/dt qvel
%u_err:     na x 1, error for actuator model
%out:       ForwardDynOutput object


%OPTIONS
do_reuse_H = 0;
do_approx_C = 0; %treat WMR as single rigid body when computing joint space bias force

if mdl.use_constraints && ~isempty(mdl.wgc_fh)
    %force balance optimization parameters
    qr_tol = 1e-3;
    max_iter = 20;
    cost_tol = 1;
    dcost_tol = cost_tol/10;
end

%use persistent for variables that only need to be computed once
%abbreviate some WmrModel properties by copying them to independent vars
persistent H inv_H Ist

%get from WmrModel object
ns = length(state0);
nv = mdl.nf+5; %number of elements in qvel

[~,~,isjd] = stateIs(ns);

%contacts
incontact = [contacts.incontact];
np = length(incontact);
npic = sum(incontact);

%to compute wheel/sprocket radius * angular rate
if mdl.nw > 0
    wsframeinds = mdl.wheelframeinds; %wheel/sprocket frame indices
    whichwt = 1:mdl.nw; %which wheel/track does contact point belong to
elseif mdl.nt > 0
    wsframeinds = mdl.sprocketframeinds;
    whichwt = uncatIndexArray([contacts.np]);
end
radii = [mdl.frames(wsframeinds).rad];

%determine which DOF are actuated from ucmd vs. WmrModel

%'actuated' can mean 'sensed' in estimation context
vis_act = ~isnan(u_cmd)'; %1 x nv, which elements of joint space velocity (qvel) are actuated
u_cmd = u_cmd(vis_act);
na = sum(vis_act); %number of actuated DOF


%init outputs that may not be set later
u_err = zeros(na,1);
out = ForwardDynOutput();
out.vc0 = NaN(3,np);
out.modelf_contact = zeros(3,np);
out.modelf_act = zeros(na,1);
out.iter = 0;
out.cost = 0;


%spatial transforms
Xup = invHTsToPluckers(HT_parent);
%Xup(:,:,i) transforms spatial *motion* vector from frame parent(i) to i coords
%Xup(:,:,i)' transforms spatial *force* vector from frame i to parent(i) coords!



%joint space inertia
if ~do_reuse_H || isempty(H)

    %Ist(:,:,i) is the 6x6 spatial inertia of the subtree rooted at frame i, in the coordinates of frame i
    Ist = subtreeInertias(mdl,Xup);
    H = jointSpaceInertia(mdl,Xup,Ist);
end

if ~do_reuse_H || isempty(inv_H)
    %invert H
    %matrix inversion is slow & innaccurate in general, but inv(H) is required many times so precomputing saves time?
    %faster than inv(H)?

    U = chol(H);
    invU = U\eye(size(U)); %backslash operator recognizes upper triangular systems
    inv_H = invU*invU';
end


%joint space bias force
dz0 = [contacts.dz];
if ~do_approx_C
    C = jointSpaceBiasForce(mdl,Xup,qvel0);
else
    %approximate, all inertia at body frame
    %spatial acc of gravity in body coords
    g = zeros(6,1);
    g(4:6) = Xup(1:3,3,1)*-mdl.grav; %world z axis in body coords * grav acc

    C = zeros(nv,1);
    C(1:6) = -Ist(:,:,1)*g + crf(qvel0(1:6))*Ist(:,:,1)*qvel0(1:6);
end


tau = zeros(nv,1);

%count constraints
ncc = 0;
if npic > 0
    if mdl.nw > 0
        Acontact = wheelJacobians(mdl, HT_world, contacts);
    elseif mdl.nt > 0
        Acontact = trackJacobians(mdl, HT_world, contacts);
    end
    
    vc0_incontact = reshape(Acontact*qvel0,3,npic);
    out.vc0(:,incontact) = vc0_incontact;
    
    if mdl.use_constraints
        ncc=3*npic; %number of contact constraints
    else
        Rw0 = radii(whichwt) .* qvel0(wsframeinds(whichwt)+5)';
        fc_ = feval(mdl.wgc_fh,mdl.wgc_p,incontact,vc0_incontact,Rw0(incontact),dz0(incontact));
        tau = tau + Acontact'*fc_(:);
        
        out.modelf_contact(:,incontact) = fc_;
    end
    
end


if mdl.use_constraints
    nac = na; %number of actuator constraints
else
    %TODO, this is unstable!
    nac=0;
    [out.modelf_act,u_err] = feval(mdl.act_fh, mdl.act_p, u_cmd, qvel0(vis_act), u_interr);
    tau(vis_act) = tau(vis_act) + out.modelf_act;
end

njc=0;
jd0 = state0(isjd);
if ~isempty(mdl.hjc_fh)
    [c,Jc,out.modelf_joint] = feval(mdl.hjc_fh, mdl, jd0, qvel0(7:end));
    
    if mdl.use_constraints
        njc=length(out.modelf_joint); %number of holonomic joint constraints
        Ajoint = [zeros(njc,6), Jc]; %Jacobian wrt qvel
    else
        njc=0;
        tau(7:end) = tau(7:end) + Jc'*out.modelf_joint;
    end
    
end

nc = ncc+nac+njc; %total number of constraints

if nc==0
    %no constraints shortcut
    qacc = H\(tau-C);
    return
end


%minimize the norm of H*qacc-(tau-C) or equivalently:
%min 0.5*qacc'*H*qacc - (tau-C)'*qacc
%subject to:
%A*qacc=b (wheel-ground contact constraints, actuator constraints, holonomic joint constraints)
%solve using Lagrange multipliers:
%[H A'; A 0]*[qacc; lambda] = [tau-C; b];
%First, solve for lambda (i.e. constraint forces/torques), then solve for qacc.
%must handle overconstraint, force balance optimization to solve for b

%allocate
A = zeros(nc,nv); %A*qacc=b
x = NaN(nc,1); %input vector for force balance optimization

%wheel constraints
cis.contact = false(1,nc); %is wheel constraint

if ncc > 0
    cis.contact(1:ncc) = true; 

    A(cis.contact,:) = Acontact;

    tmp = vc0_incontact;
    tmp(3,:) = 0; %init vz to zero, or never reaches steady state

    x(cis.contact) = tmp(:);
end

%actuator constraints
cis.act = false(1,nc); %is actuator constraint
if nac > 0
    cis.act((ncc+1):(ncc+nac)) = true; 
    Aact = eye(nv); Aact = Aact(vis_act,:);
    A(cis.act,:) = Aact;
    x(cis.act)=qvel0(vis_act);
%     x(cis.act)=ucmd; %DEBUGGING
end

%holonomic joint constraints
cis.joint = false(1,nc); %is joint constraint
if njc > 0;
    cis.joint((ncc+nac+1):end)=true;
    A(cis.joint,:) = Ajoint;
    x(cis.joint) = Ajoint*qvel0;
end


if isempty(mdl.wgc_fh)
    useErpCfm()
    
    return
end


%Get a maximal linearly independent subset of constraints (rows of A)
cis.ind = true(1,nc);


if (nc-nac) > 0
    %actuator constraints are automatically included
    A_ = A(~cis.act,:);
    A_(:,vis_act) = [];

    [~, R, E] = qr(A_',0); %A_' is (nv-na) x (nc-nac)
    diagR = abs(diag(R));
    is_ind_ = false(nc-nac,1); %logical, rows of A_ in linearly independent subset
    is_ind_(E(diagR >= qr_tol*diagR(1))) = true; %length(diagR) < length(E) okay
    
    cis.ind(~cis.act) = is_ind_;

end
cis.inpt = cis.ind;

if isempty(mdl.act_fh)
    cis.inpt(cis.act)=false; 
end %ideal actuators

x=x(cis.inpt);
ninpt=length(x);

Aind = A(cis.ind,:);

%precompute vars required for cost function

%compute Hl,fl such that Hl*lambda = fl;
Hl = Aind*inv_H*Aind';

%invert Hl
% inv_Hl = inv(Hl);
%Hl is symmetric, this is faster (?)
Ul = chol(Hl);
invUl = Ul\eye(size(Ul));
inv_Hl = invUl*invUl';
        
fl_c = - Aind*(inv_H*(tau-C)); %constant part of fl

%initialize variables in nested functions to avoid dynamic assignment
%computed in calcCost, needed in calcGradient
mf = zeros(nc,1); %model force
dmf = struct; %derivative of model force
err = []; %force-balance error
dmf_dx = []; %DEBUGGING

%TODO, move this block to force_balance function?
fh_calc_cost = @calcCost;
fh_calc_gradient = @calcGradient;

dt_rem = dt; %remaining step size

%DEBUGGING
backup = struct;
backup.dz0 = dz0;
backup.jd0 = jd0;
backup.qvel0 = qvel0;
backup.dt = dt;

iter = 0;

while 1
    
    cost = calcCost(x);

    iter_ = 0;
%     fprintf('iter: %d, cost: %f\n', iter_, cost) %DEBUGGING

%     if 1 %DEBUGGING, at least one Newton's method step
    if cost > cost_tol
        %optimization required

        [dqvel_dx, dtaulambda_dx, dwgcin_dx] = precompute();
        grad = calcGradient;

        cost_prev = cost;

        for iter_ = 1:max_iter

%             calcGradientNum(x)

            Hess = 2*(derr_dx'*derr_dx);
            p = -Hess\grad'; %direction
            
            [x,cost,grad] = linesearch(x,cost,grad,p,1,fh_calc_cost,fh_calc_gradient);
            
%             fprintf('iter: %d, cost: %f\n', iter_, cost) %DEBUGGING

            if cost < cost_tol %global minimum
                break
            end

            if cost_prev-cost < dcost_tol %local minimum
                if 0 %partial step
                                
                    dt_prev = dt; %backup
                    dt = dt/2;

%                     [dqvel_dx, dtaulambda_dx, dwgcin_dx] = precompute();
                    changeDt(dt_prev); %faster!

                    cost=calcCost(x);
                    grad=calcGradient;
                    cost_prev=Inf;
                    fprintf('decreased step size to %.4f\n',dt) %DEBUGGING
                else
                    break
                end
            else
                cost_prev = cost;
            end
        end
    end
    
    dt_rem = dt_rem - dt;
    iter = iter + iter_;
    
    if dt_rem == 0
        break
    end
    
    %take a partial step, integrate
    if nwc > 0, dz0(incontact) = dz_incontact; end
    if njc > 0, jd0 = jd; end
    qvel0 = qvel;
    dt = dt_rem;
    
    fprintf('remaining step size: %.4f\n',dt) %DEBUGGING
    
end

qacc = (qvel - backup.qvel0)/backup.dt;

%outputs
out.iter = iter;
out.cost = cost;

if ncc > 0
    out.modelf_contact(:,incontact) = reshape(mf(cis.contact),3,npic);
end

if nac > 0 && ~isempty(mdl.act_fh)
    out.modelf_act = mf(cis.act);
end

if njc > 0
    out.modelf_joint = mf(cis.joint);
end

%DEBUGGING
out.vc = NaN(3,np);
if ncc > 0
    out.vc(:,incontact) = vc_incontact;
end

    function [dqvel_dx, dtaulambda_dx, dwgcin_dx] = precompute()
        
        %precompute vars required for gradient function
        %init output
        dwgcin_dx = [];
        
        dlambda_dx = inv_Hl(:,cis.inpt(cis.ind))/dt; %dlambda/dx
        dqvel_dx = (inv_H*Aind'*dt)*dlambda_dx; %dqvel/dlambda*dlambda/dx
        
        if ncc > 0
            %d(wheel-ground contact model inputs)/dx
            dwgcin_dx = zeros(5,ninpt,npic);
            
            %TODO, simplify?
            picno = 0;
            for pno = find(incontact)
                wtno = whichwt(pno); %wheel/sprocket number
                picno = picno+1; %point in contact number
                rows = (1:3) + (picno-1)*3;
                dwgcin_dx(1:3,:,picno) = Acontact(rows,:)*dqvel_dx;
                dwgcin_dx(4,:,picno) = radii(wtno)*dqvel_dx(wsframeinds(wtno)+5,:);
                dwgcin_dx(5,:,picno) = dwgcin_dx(3,:,picno)*dt;
            end
        end
        dtaulambda_dx = Aind'*dlambda_dx;
        
    end

    function changeDt(dt_prev)
        %faster than calling precompute()
        %dt_prev is previous step size
        %dqvel_dx unchanged

        s = dt/dt_prev;
        dtaulambda_dx = dtaulambda_dx/s;
        dwgcin_dx(5,:,:) = dwgcin_dx(5,:,:)*s;
    end

    function cost = calcCost(x_)
        xb=NaN(nc,1); %buffered
        xb(cis.inpt)=x_;
        
        b=NaN(nc,1);
        if ncc > 0
            b(cis.contact) = xb(cis.contact) - vc0_incontact(:);
        end
        if nac > 0
            if isempty(mdl.act_fh)
                b(cis.act) = u_cmd - qvel0(vis_act); %ideal actuators
            else
                b(cis.act) = xb(cis.act) - qvel0(vis_act);
            end
        end
        if njc > 0
            b(cis.joint) = xb(cis.joint) - Ajoint*qvel0; 
        end
        
        fl = b(cis.ind)/dt + fl_c;
        
        %constraint forces
        lambda = inv_Hl*fl;
        qacc = inv_H*(tau-C+Aind'*lambda); 
        qvel = qvel0+qacc*dt; %new qvel, time i+1
        
        %wheel-ground contact model
        if ncc > 0
            vc_incontact = reshape(Acontact*qvel,3,npic); %velocity of contact points at time i+1, for wheels in contact
            Rw = radii(whichwt) .* qvel(wsframeinds(whichwt)+5)';
            dz_incontact = dz0(incontact) + vc_incontact(3,:)*dt; %symplectic Euler
            
            [mf(cis.contact),dmf.wheel_dwgcin] = feval(mdl.wgc_fh,mdl.wgc_p,incontact,vc_incontact,Rw(incontact),dz_incontact);
        end

        %actuator model
        if nac > 0 && ~isempty(mdl.act_fh);
            [mf(cis.act), u_err, dmf.act_du] = feval(mdl.act_fh, mdl.act_p, u_cmd, qvel(vis_act), u_interr);
        end
        
        %holonomic joint constraint model
        if njc > 0
            jd = jd0 + qvel(7:end)*dt; %symplectic Euler
            [~, ~, mf(cis.joint), dmf.joint_djd, dmf.joint_djr] = feval(mdl.hjc_fh,mdl,jd,qvel(7:end));
        end
        
        tau_lambda = Aind'*lambda;
        tau_mf = A'*mf;

        err = tau_lambda - tau_mf; %error
        
        if nac > 0 && isempty(mdl.act_fh) %ideal actuators
            err(vis_act) = 0;
        end
        
        cost = err'*err;

    end

    function grad_ = calcGradient()

        dmf_dx = zeros(nc,ninpt);
        
        if ncc > 0
            dmfcontact_dx = zeros(ncc,ninpt);
            for picno = 1:npic
                rows = (1:3) + (picno-1)*3;
                dmfcontact_dx(rows,:) = dmf.wheel_dwgcin(:,:,picno)*dwgcin_dx(:,:,picno); %dfc/dwgcin*dwgcin/dx
            end
            dmf_dx(cis.contact,:) = dmfcontact_dx;
        end
        if nac > 0 && ~isempty(mdl.act_fh)
            dmfact_dx = dmf.act_du*dqvel_dx(vis_act,:);
            dmf_dx(cis.act,:) = dmfact_dx;
        end
        if njc > 0
            djd_dx = dqvel_dx(7:end,:)*dt;
            dmfjoint_dx = dmf.joint_djd*(djd_dx) + dmf.joint_djr*dqvel_dx(7:end,:);
            dmf_dx(cis.joint,:) = dmfjoint_dx;
        end
        dtaumf_dx = A'*dmf_dx;
        
        derr_dx = dtaulambda_dx - dtaumf_dx; %required to compute Hessian
        
        if nac > 0 && isempty(mdl.act_fh) %ideal actuators
            derr_dx(vis_act,:) = 0;
        end
        
        grad_ = 2*err'*derr_dx;
        
    end

    function calcGradientNum(x)
        %DEBUGGING, compute the gradient numerically
        
        x_backup = x;
        cost_backup = cost;
%         err_backup = err;
%         mf_backup = mf;
        grad_num = zeros(1,ninpt);
%         derr_dx_num = zeros(length(err),ninpt);
%         dmf_dx_num = zeros(nc,ninpt);
        myeps = 1e-8; %how to choose?
        for k = 1:ninpt
            x = x_backup;
            x(k) = x(k) + myeps;
            cost = calcCost(x);
            grad_num(k) = (cost-cost_backup)/myeps;
%             derr_dx_num(:,k) = (err-err_backup)/myeps;
%             dmf_dx_num(:,k) = (mf-mf_backup)/myeps;
        end
        x = x_backup;
        cost = calcCost(x); %to restore cost and other vars
        
        %print normalized error
        fprintf('norm. grad error: %s\n',num2str((grad-grad_num)/norm(grad_num)*100,'%.4f '))
        
    end

    function useErpCfm()
        %use erp, cfm like Open Dynamics Engine
        %TODO, check this
        %TODO, limits on contact forces. (normal force > 0, tangential force < mu*normal force)
        
        b = NaN(nc,1);
        cfm_diag = zeros(nc,1);

        if ncc > 0
            vc = zeros(3,np);
            vc(3,:) = -mdl.erp_z(whichwt) .* dz0/dt;
            vc_incontact = vc(:,incontact);
            temp = vc_incontact - vc0_incontact;

            b(cis.contact) = temp(:);

            temp = [mdl.fds_x(whichwt); mdl.fds_y(whichwt); mdl.cfm_z(whichwt)];
            temp = temp(:,incontact);
            cfm_diag(cis.contact) = temp(:);
        end
        if nac > 0
            b(cis.act) = u_cmd - qvel0(vis_act);
        end
        if njc > 0
            b(cis.joint) = -mdl.erp_j .* c/dt - Ajoint*qvel0;
            cfm_diag(cis.joint) = mdl.cfm_j;
        end
        b = b/dt;


        %Hl*lambda = fl
        Cfm = diag(cfm_diag/dt); %why divide by dt? because solving for acceleration (not velocity)
        Hl = A*inv_H*A' + Cfm;
        fl = b - A*inv_H*(tau-C);

        %TODO, constraints on lambda, solve QP
        %-normal force >= 0
        %-lon,lat force <= lim

        lambda = Hl\fl;

        qacc = inv_H*(tau-C+A'*lambda);

        %outputs
        if ncc > 0
            out.modelf_contact(:,incontact) = reshape(lambda(cis.contact),3,npic);
        end

        if nac > 0
            out.modelf_act = lambda(cis.act);
        end

        if njc > 0
            out.modelf_joint = lambda(cis.joint);
        end

        out.vc = NaN(3,np);
        if ncc > 0
            %account for cfm
            vc_incontact = vc_incontact - reshape(cfm_diag(cis.contact).*lambda(cis.contact),3,npic);
            out.vc(:,incontact) = vc_incontact;
        end
    end

end


