%test_kinematics.m

%%

%test_wheelJacobians
clear_all
close all
clc

[mdl, state] = zoe();

surfs = {};
surfs{end+1} = flat();

%convert state to homogeneous transforms
[~,HT_world] = stateToHT(mdl,state);

%contact geometry
contacts(1:mdl.nw) = WheelContactGeom();
contacts = updateModelContactGeom(mdl, surfs, HT_world, mdl.min_npic, contacts);
    
A = wheelJacobians(mdl,HT_world,contacts)

if 1
    % time it
    irf = 1e-3;
    n = 1e6;
    
    tic
    for i=1:(n*irf)
        A = wheelJacobians(mdl,HT_world,contacts);
    end
    t=toc;

    fprintf('iterations: %e\n',n)
    t_=t/irf*1e3;
    fprintf('Elapsed time (ms): %f\n',t_)
end

%%
% test_trackJacobians
clear_all
close all
clc

[mdl, state] = talon();

surfs = {};
surfs{end+1} = flat();

%convert state to homogeneous transforms
[~,HT_world] = stateToHT(mdl,state);

%contact geometry
contacts = initTrackContactGeom(mdl.frames(mdl.sprocketframeinds));
contacts = updateModelContactGeom(mdl, surfs, HT_world, mdl.min_npic, contacts);

A = trackJacobians(mdl,HT_world,contacts)

if 1
    % time it
    irf = 1e-3;
    n = 1e6;
    
    tic
    for i=1:(n*irf)
        A = trackJacobians(mdl,HT_world,contacts);
    end
    t=toc;

    fprintf('iterations: %e\n',n)
    t_=t/irf*1e3;
    fprintf('Elapsed time (ms): %f\n',t_)
end



%%
%test_forwardVelKin

clear_all
close all
clc

%make WmrModel object
[mdl, state] = zoe();
% [mdl, state] = rocky();
% [mdl, state] = talon();

nv = mdl.nf+5;

mdl.bsm_fh = [];

%terrain
surfs = {};
surfs{end+1} = flat();

%init contact geometry
if mdl.nw > 0
    contacts(1:mdl.nw) = WheelContactGeom();
elseif mdl.nt > 0
    contacts = initTrackContactGeom(mdl.frames(mdl.sprocketframeinds));
end

%***copied from odeKin
%convert state to homogeneous transforms
[~,HT_world] = stateToHT(mdl,state);

%update contact geometry
contacts = updateModelContactGeom(mdl, surfs, HT_world, mdl.min_npic, contacts);

%get control inputs
[u,qvel_cmd] = feval(mdl.controller_fh,mdl,0,state)
vis_act = [];

%***end copied section

[qvel,vc] = forwardVelKin(mdl,state,u,vis_act,HT_world,contacts)


if 1
    % time it
    irf = 1e-2;
    n = 1e5;
    
    tic
    for i=1:(n*irf)
        qvel = forwardVelKin(mdl,state,u,vis_act,HT_world,contacts);
    end
    t=toc;


    fprintf('iterations: %e\n',n)
    t_=t/irf*1e3;
    fprintf('Elapsed time (ms): %f\n',t_)
end




%%
%test_initTerrainContact

clear_all
close all
clc

do_anim = 1;

%make WmrModel object
model_fh = @zoe;
% model_fh = @rocky;
% model_fh = @talon;

%make WmrModel object
if do_anim
    [mdl, state, qvel, anim] = feval(model_fh);
else
    [mdl, state, qvel] = feval(model_fh);
end

%terrain
surfs = {};

% surfs{end+1} = flat();
% surfs{end+1} = ramp();

setseed(123);
surfs{end+1} = randomgrid();
printGridSurf(surfs{end},[resourcedir() 'gridsurfdata.txt'])

%init contact geometry
if mdl.nw > 0
    contacts(1:mdl.nw) = WheelContactGeom();
elseif mdl.nt > 0
    contacts = initTrackContactGeom(mdl.frames(mdl.sprocketframeinds));
end

state_before = state
state = initTerrainContact(mdl,surfs,contacts,state_before);

state

if 1
    %time it
    irf = 1e-2;
    n = 1e4;

%     profile on
    tic
    for i=1:(n*irf)
        dummy = initTerrainContact(mdl,surfs,contacts,state_before);
    end
    t=toc;
%     profile off
%     profile viewer

    fprintf('iterations: %e\n',n)
    t_=t/irf*1e3;
    fprintf('Elapsed time (ms): %f\n',t_)
end



if do_anim
    %%

    drawSurfaces(surfs);
    
    %copied from odeKin
    %convert state to homogeneous transforms
    [HT_parent, HT_world] = stateToHT(mdl,state);
    
    %update contact geometry
    contacts = updateModelContactGeom(mdl, surfs, HT_world, mdl.min_npic, contacts);
    
    updateAnimation(anim, HT_parent, contacts);
    centerAxes(anim, HT_parent(:,:,1));
    
    zoom(1.5)
    resizeFig(2,2)
    view(30,30)

end



    
