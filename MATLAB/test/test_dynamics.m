%test_dynamics.m

%%

%test_subtreeInertias

clear_all
close all
clc

[mdl, state] = zoe();
% [mdl, state] = rocky();

HT_parent = stateToHT(mdl,state);

% Is = mdl.Is
Xup = invHTsToPluckers(HT_parent);
Is_subt = subtreeInertias(mdl,Xup)

if 1
    % time it
    irf = 1e-2;
    n = 1e6;
    tic
    for i=1:(n*irf)
        Is_subt = subtreeInertias(mdl,Xup);
    end
    t=toc;

    fprintf('iterations: %e\n',n)
    t_=t/irf*1e3;
    fprintf('Elapsed time (ms): %f\n',t_)
end

%%
%test_jointSpaceInertia

clear_all
close all
clc

[mdl, state] = zoe();
% [mdl, state] = rocky();

[HT_parent,HT_world] = stateToHT(mdl,state);

Xup = invHTsToPluckers(HT_parent);
Is_subt = subtreeInertias(mdl,Xup);

H = jointSpaceInertia(mdl,Xup,Is_subt)

if 1
    % time it
    irf = 1e-2;
    n = 1e6;
    tic
    for i=1:(n*irf)
        H = jointSpaceInertia(mdl,Xup,Is_subt);
    end
    t=toc;

    fprintf('iterations: %e\n',n)
    t_=t/irf*1e3;
    fprintf('Elapsed time (ms): %f\n',t_)
end

%%
%test_jointSpaceBiasForce

clear_all
close all
clc

[mdl, state, qvel] = zoe();
% [mdl, state, qvel] = rocky();

[HT_parent,HT_world] = stateToHT(mdl,state);

Xup = invHTsToPluckers(HT_parent);

qvel = (1:length(qvel))' %DEBUGGING

C = jointSpaceBiasForce(mdl, Xup, qvel);

disp('C = ')
disp(num2str(C,'%.4f'))

if 1
    % time it
    irf = 1e-3;
    n = 1e6;
    tic
    for i=1:(n*irf)
        C = jointSpaceBiasForce(mdl, Xup, qvel);
    end
    t=toc;

    fprintf('iterations: %e\n',n)
    t_=t/irf*1e3;
    fprintf('Elapsed time (ms): %f\n',t_)
end


%%
%test_forwardDyn

clear_all
close all
clc

[mdl, state, qvel] = zoe();
% [mdl, state, qvel] = rocky();
% [mdl, state, qvel] = talon();

nv = mdl.nf + 5;

% qvel(:) = 0; 
mdl.act_fh = []; %ideal actuators

surfs = {};
surfs{end+1} = flat();

%init contact geometry
if mdl.nw > 0
    contacts(1:mdl.nw) = WheelContactGeom();
elseif mdl.nt > 0
    contacts = initTrackContactGeom(mdl.frames(mdl.sprocketframeinds));
end

% state = initTerrainContact(mdl,surfs,contacts,state);

%***copied from odeDyn
%convert state to homogeneous transforms
[HT_parent, HT_world] = stateToHT(mdl,state);

%update contact geometry
contacts = updateModelContactGeom(mdl, surfs, HT_world, 0, contacts);

%get control inputs
u = ControllerIO();

u.cmd = feval(mdl.controller_fh,mdl,0,state);

u.vis_act = false(1,nv); 
u.vis_act(mdl.actframeinds+5) = true;

%***end copied section

u.interr = zeros(mdl.na,1);
dt = .04;

state
qvel
% u_
qacc = forwardDyn(mdl, state, qvel, u, HT_parent, HT_world, contacts, dt)


if 1
    % time it
    irf = 1e-3;
    n = 1e5;
    
    tic
    for i=1:(n*irf)
        qacc = forwardDyn(mdl, state, qvel, u, HT_parent, HT_world, contacts, dt);
    end
    t=toc;


    fprintf('iterations: %e\n',n)
    t_=t/irf*1e3;
    fprintf('Elapsed time (ms): %f\n',t_)
end