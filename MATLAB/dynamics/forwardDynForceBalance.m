function [qacc, u, out] = forwardDynForceBalance(mdl, state0, qvel0, u, contacts, dt, H, C, A, cis)
%Compute *constrained* forward dynamics using force-balance optimization

%OPTIONS
% qr_tol = 1e-3;
max_iter = 20;
err_tol = 1e-4;
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

jd0 = state0(isjd);

%END PRE-PROCESS

tau = zeros(nv,1);

%allocate

%if ideal actuators
if nac > 0 && isempty(mdl.act_fh) 
    %ideal actuators
    ninpt = nv - nac;
else
    ninpt = nv; %qacc
end
x = zeros(ninpt,1); %input vector for force-balance optimization


%initialize variables in nested functions to avoid dynamic assignment
%computed in calcCost, needed in calcGradient
mf = zeros(nc,1); %model force
dmf = struct; %derivative of model force
err = []; %force-balance error
% dmf_dqacc = []; %DEBUGGING


fh_calc_cost = @calcCost;
fh_calc_gradient = @calcGradient;
    
cost = calcCost(x);

iter = 0;

% fprintf('iter: %d, max(abs(err)): %f, cost: %f\n', iter, max(abs(err)), cost) %DEBUGGING

if max(abs(err)) > err_tol
    %optimization required

    [dwgcin_dqacc] = precomputeGrad();
    grad = calcGradient;

    cost_prev = cost;

    for iter = 1:max_iter

%         calcGradientNum(x)

        Hess = 2*(derr_dx'*derr_dx);
        
        [Hess_U, Hess_p] = chol(Hess); %Hess_U is upper triangular
        if Hess_p == 0
            %Hess is positive definite, solve using Cholesky decomposition (faster)
            p = -Hess_U\(Hess_U'\grad');
%             p = -Hess\grad'; %DEBUGGING
        else %else use pseudo-inverse, TODO CHECK THIS
            p = -pinv(Hess)*grad';
        end

        
%         p = -Hess\grad'; %direction

        [x,cost,grad] = linesearch(x,cost,grad,p,1,fh_calc_cost,fh_calc_gradient);

%         fprintf('iter: %d, max(abs(err)): %f, cost: %f\n', iter, max(abs(err)), cost) %DEBUGGING

        if max(abs(err)) <= err_tol
            break
        end
        
%         if cost <= cost_tol %global minimum
%             break
%         end

        if cost_prev-cost <= dcost_tol %local minimum
            break
        end
        cost_prev = cost;
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

    function [dwgcin_dqacc] = precomputeGrad()
        
        %precompute vars required for gradient function, 
        %derivatives with respect to inputs x
        
        %init output
        %TODO, eliminate or simplify this variable
        dwgcin_dqacc = [];
        
        if ncc > 0
            dwgcin_dqacc = zeros(5,nv,npic);
            picno = 0;
            for pno = find(incontact)
                wtno = whichwt(pno); %wheel/track number
                picno = picno+1; %point in contact number
                rows = (1:3) + (picno-1)*3;
                dwgcin_dqacc(1:3,:,picno) = Acontact(rows,:)*dt; %vx vy vz
                dwgcin_dqacc(4,wsframeinds(wtno)+5,picno) = radii(wtno)*dt; %R*omega
                dwgcin_dqacc(5,:,picno) = dwgcin_dqacc(3,:,picno)*dt; %dz
            end
        end
        
    end

    function cost = calcCost(x_)
        
        if nac > 0 && isempty(mdl.act_fh)
            %ideal actuators
            qacc = zeros(nv,1);
            qacc(~vis_act) = x_;
            qacc(vis_act) = (u.cmd - qvel0(vis_act))/dt;
        else
            qacc = x_;
        end
        qvel = qvel0 + qacc*dt;
        tau_lambda = -(tau - C) + H*qacc;
        
        %compute model forces
        %wheel-ground contact model
        if ncc > 0
            vc_incontact = reshape(Acontact*qvel,3,npic); %velocity of contact points at time i+1, for wheels in contact
            Rw = radii(whichwt) .* qvel(wsframeinds(whichwt)+5)';
            dz_incontact = dz0(incontact) + vc_incontact(3,:)*dt; %symplectic Euler
            
            [mf(cis.contact),dmf.contact_dwgcin] = feval(mdl.wgc_fh,mdl.wgc_p,incontact,vc_incontact,Rw(incontact),dz_incontact);
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
        
        tau_mf = A'*mf;          
        
        err = tau_lambda - tau_mf; %error
        
        if nac > 0 && isempty(mdl.act_fh) %ideal actuators
            err(vis_act) = 0;
        end
        
        cost = err'*err;
        
%         blah = 1;

    end

    function grad_ = calcGradient()
        dmf_dqacc = zeros(nc,nv);

        if ncc > 0
            off = find(cis.contact,1)-1; %offset, based on index of first contact constraint
            for picno = 1:npic
                rows = (1:3) + (picno-1)*3 + off;
                dmf_dqacc(rows,:) = dmf.contact_dwgcin(:,:,picno)*dwgcin_dqacc(:,:,picno);
            end
        end
        if nac > 0 && ~isempty(mdl.act_fh)
            dmf_dqacc(cis.act,vis_act) = dmf.act_du*dt;
        end
        if njc > 0
            dmf_dqacc(cis.joint,7:end) = (dmf.joint_djd * dt^2) + (dmf.joint_djr * dt);
        end

        if nac > 0 && isempty(mdl.act_fh)
            %ideal actuators
            dtaulambda_dx = H(:,~vis_act);
            dtaumf_dx = A'*dmf_dqacc(:,~vis_act);
        else
            dtaulambda_dx = H;
            dtaumf_dx = A'*dmf_dqacc;
        end

        derr_dx = dtaulambda_dx - dtaumf_dx; %required to compute Hessian
        
        if nac > 0 && isempty(mdl.act_fh) %ideal actuators
            derr_dx(vis_act,:) = 0;
        end
        
        grad_ = 2*err'*derr_dx;
        
%         blah = 1;
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
        
        myeps = 1e-6; %how to choose?
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
        
%         blah = 1;
        
    end

   
end

