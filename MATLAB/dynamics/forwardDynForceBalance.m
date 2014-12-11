function [qacc, u, out] = forwardDynForceBalance(mdl, state0, qvel0, u, contacts, dt, inv_H, C, A, cis)
%Compute *constrained* forward dynamics using force-balance optimization

%OPTIONS
qr_tol = 1e-3;
max_iter = 20;
cost_tol = 1;
dcost_tol = cost_tol/10;

%PRE-PROCESS

%get from WmrModel object
ns = length(state0);
nv = mdl.nf+5; %number of elements in qvel
[~,~,isjd] = stateIs(ns);

%to compute wheel/sprocket radius * angular rate
if mdl.nw > 0
    wsframeinds = mdl.wheelframeinds; %wheel/sprocket frame indices
    whichwt = 1:mdl.nw; %which wheel/track does contact point belong to
elseif mdl.nt > 0
    wsframeinds = mdl.sprocketframeinds;
    whichwt = uncatIndexArray([contacts.np]);
end
radii = [mdl.frames(wsframeinds).rad];

%contacts
incontact = [contacts.incontact];
np = length(incontact);
npic = sum(incontact);
dz0 = [contacts.dz];

%actuators
vis_act = u.vis_act; %abbreviate
na = sum(vis_act); %number of actuated DOF

%init output
u.err = zeros(na,1); %in case ideal actuators

out = ForwardDynOutput();
out.vc0 = NaN(3,np);
out.modelf_contact = zeros(3,np);
out.modelf_act = zeros(na,1); %in case ideal actuators

%
nc = size(A,1);
ncc = sum(cis.contact);
nac = sum(cis.act);
njc = sum(cis.joint);

Acontact = A(cis.contact,:);
Ajoint = A(cis.joint,:);

vc0_incontact = reshape(A(cis.contact,:)*qvel0, 3, npic);

jd0 = state0(isjd);

%END PRE-PROCESS

tau = zeros(nv,1);

%allocate
x = NaN(nc,1); %input vector for force-balance optimization

%wheel constraints
if ncc > 0
    tmp = vc0_incontact;
    tmp(3,:) = 0; %init vz to zero, or never reaches steady state

    x(cis.contact) = tmp(:);
end

%actuator constraints
if nac > 0
    x(cis.act) = qvel0(vis_act);
%     x(cis.act) = u.cmd; %DEBUGGING
end

%holonomic joint constraints
if njc > 0;
    x(cis.joint) = Ajoint*qvel0;
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

x = x(cis.inpt);
ninpt = length(x);

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
% dmf_dx = []; %DEBUGGING


fh_calc_cost = @calcCost;
fh_calc_gradient = @calcGradient;
    
cost = calcCost(x);

iter = 0;

% fprintf('iter: %d, cost: %f\n', iter, cost) %DEBUGGING

if cost > cost_tol
    %optimization required

    [dqvel_dx, dwgcin_dx, dtaulambda_dx] = precomputeGrad();
    grad = calcGradient;

    cost_prev = cost;

    for iter = 1:max_iter

%         calcGradientNum(x)

        Hess = 2*(derr_dx'*derr_dx);
        p = -Hess\grad'; %direction

        [x,cost,grad] = linesearch(x,cost,grad,p,1,fh_calc_cost,fh_calc_gradient);

%         fprintf('iter: %d, cost: %f\n', iter, cost) %DEBUGGING

        if cost < cost_tol %global minimum
            break
        end

        if cost_prev-cost < dcost_tol %local minimum
            break
        else
            cost_prev = cost;
        end
    end
end
    

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

    function [dqvel_dx, dwgcin_dx, dtaulambda_dx] = precomputeGrad()
        
        %precompute vars required for gradient function, 
        %derivatives with respect to inputs x
        
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

    function cost = calcCost(x_)
        xb=NaN(nc,1); %buffered
        xb(cis.inpt)=x_;
        
        b=NaN(nc,1);
        if ncc > 0
            b(cis.contact) = xb(cis.contact) - vc0_incontact(:);
        end
        if nac > 0
            if isempty(mdl.act_fh)
                b(cis.act) = u.cmd - qvel0(vis_act); %ideal actuators
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
            [mf(cis.act), u.err, dmf.act_du] = feval(mdl.act_fh, mdl.act_p, u.cmd, qvel(vis_act), u.interr);
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

   
end

