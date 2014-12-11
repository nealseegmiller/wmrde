%test_state.m



%%
% test_stateToHT()
clear_all
close all
clc

% [mdl, state] = rocky();
[mdl, state] = zoe();

state
[HT_parent,HT_world] = stateToHT(mdl,state)

if 1
    % time it
    irf = 1e-2;
    n = 1e6;
    tic
    for i=1:(n*irf)
        [HT_parent,HT_world] = stateToHT(mdl,state);
    end
    t=toc;

    fprintf('iterations: %e\n',n)
    t_=t/irf*1e3;
    fprintf('Elapsed time (ms): %f\n',t_)
end



%%
% test_qvelToQdot()
clear all
close all
clc

olen=3;
nf=5;
ns=nf-1+olen+3;
nv=nf-1+3+3;

euler = [10 20 30]'*pi/180;

orient = euler
% orient = eulerToQuat(euler)

R_body_to_world = eulerToRot(euler)

qvel = (1:nv)'-1

statedot = qvelToQdot(qvel,euler,R_body_to_world)

qvel_convertedback = qdotToQvel(statedot,euler,R_body_to_world)


