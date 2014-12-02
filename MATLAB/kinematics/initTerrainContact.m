function state = initTerrainContact(mdl,surfaces,contacts,state)
%initialize state such that wheels are in contact with terrain
%changes elevation, attitude, passive joints
%INPUTS
%mdl:       WmrModel object, must have .isfixed field
%surfaces:  cell array of surface objects
%state:     ns x 1 vector
%contacts:  1 x nw array of WheelContactGeom objects, or
%           1 x nt array of TrackContactGeom objects
%cost_tol:  cost tolerance, omit or empty to use default
%TODO, cost_tol input, can be Inf
%OUTPUTS
%state:     (modified)



max_iter = 15;
cost_tol = 1e-6;
dcost_tol = cost_tol/10;

ns = length(state);
[isorient,ispos,isjd]=stateIs(ns);

%convert quaternion to Euler angles if necessary, so that attitude can be changed independent of yaw.
if sum(isorient) > 3
    use_quat = true; %flag to change back on exit
    
    % back up
    state0 = state; 
    isorient0 = isorient;
    
    ns = ns-1;
    state = zeros(ns,1);
    [isorient,ispos,isjd]=stateIs(ns);
    
    state(isorient) = quatToEuler(state0(isorient0));
    state(~isorient) = state0(~isorient0);
else
    use_quat = false;
end


%contacts
if mdl.nw > 0
    whichwt = 1:mdl.nw; %which wheel/track does contact point belong to
elseif mdl.nt > 0
    whichwt = uncatIndexArray([contacts.np]);
end


%to avoid dynamic assignment
HT_world = [];
err = [];
njc = [];
Jc = [];
derr_dx = []; %required for Hessian

fh_calc_cost = @calcCost;
fh_calc_gradient = @calcGradient;

%specify free elements of state
isfree = false(1,ns);
isfree(isorient) = [1 1 0]; %rol, pit
isfree(ispos) = [0 0 1]; %z
isfree(isjd) = ~[mdl.frames(2:end).isfixed];

x = state(isfree); %init x, comprises the free elements of state

cost = calcCost(x);
cost_prev = Inf;

for iter = 1:max_iter

%     fprintf('iter: %d, cost: %f, [%s]\n', iter, cost, num2str([contacts.dz],'%.4f '))

    %break if close enough
    if cost < cost_tol
        break
    end

    %break if no improvement
    if abs(cost-cost_prev) < dcost_tol
        break
    end
    cost_prev = cost;

    if iter==1
        grad = calcGradient;
    end


    Hess = 2*(derr_dx'*derr_dx);

    p = -Hess\grad'; %direction

    if abs(grad*p) < dcost_tol
        %gradient is zero at minumum, stop if close enough.
        break
    end

    [x,cost,grad] = linesearch(x,cost,grad,p,1,fh_calc_cost,fh_calc_gradient);

end

%convert Euler angles back to quaternion if necessary
if use_quat
    state0 = state; % back up
    ns = ns+1;
    state = zeros(ns,1);
    
    state(isorient0) = euler_to_quat(state0(isorient));
    state(~isorient0) = state0(isorient);
end

    function cost = calcCost(x_)
        
        state(isfree) = x_;
        [~,HT_world] = stateToHT(mdl, state);
        
        %update contact geometry
        contacts = updateModelContactGeom(mdl, surfaces, HT_world, mdl.min_npic, contacts);

        err = [contacts.dz] - mdl.dz_target(whichwt);
        err = err([contacts.incontact])';
        
        if ~isempty(mdl.hjc_fh)
            %additional constraints
            [c, Jc] = feval(mdl.hjc_fh,mdl,state(isjd));
            err = [err; c];
            njc = length(c);
            Jc = [zeros(njc,6),Jc];
        else
            njc = 0;
        end

        cost = err'*err;
    end

    function grad = calcGradient
        
        if mdl.nw > 0
            A = wheelJacobians(mdl, HT_world, contacts);
        elseif mdl.nt > 0
            A = trackJacobians(mdl, HT_world, contacts);
        end

        derrdot_dqvel = A(3:3:end,:); %d vz/d qvel
        
        if njc > 0
            derrdot_dqvel = [derrdot_dqvel; Jc];
        end
        
        derr_dstate = qvelToQdot(derrdot_dqvel',state(isorient),HT_world(1:3,1:3,1))';
        derr_dx = derr_dstate(:,isfree);
        
        grad = 2*err'*derr_dx;

    end
end

