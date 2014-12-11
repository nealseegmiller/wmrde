function [ydot,out] = odeDyn(time,y,mdl,surfaces,contacts,dt,pathno)
%ordinary differential equation for kinematic simulation
%INPUTS
%time:      scalar
%y:         [state; qvel; interr]
%mdl:       WmrModel object
%surfaces:  cell array of surface objects for terrain
%contacts:  1 x nw array of WheelContactGeom objects, or
%           1 x nt array of TrackContactGeom objects
%dt:        only required if using constraints
%OUTPUTS
%ydot:      d/dt y
%out:       OdeOutput object

if nargin < 7
    pathno = [];
    if nargin < 6
        dt = [];
    end
end

%get from WmrModel object
nf = mdl.nf;
na = mdl.na;

nv = nf+5;
ns=length(y)-nv-na;

isorient = stateIs(ns);

u = ControllerIO();


%split up y vector
[state,qvel,u.interr]=odeDynSplitVec(y,nf,na);

%convert state to homogeneous transforms
[HT_parent, HT_world] = stateToHT(mdl,state);

%update contact geometry
contacts = updateModelContactGeom(mdl, surfaces, HT_world, 0, contacts);

%get control inputs
if isempty(pathno)
    u.cmd = feval(mdl.controller_fh,mdl,time,state);
else
    u.cmd = feval(mdl.controller_fh,mdl,time,state,pathno);
end

%compute acceleration
[qacc,u,fdout] = forwardDyn(mdl, state, qvel, u, HT_parent, HT_world, contacts, dt);

if ~isempty(dt)
    %Semi-implicit or symplectic Euler, conserves energy
    %http://en.wikipedia.org/wiki/Semi-implicit_Euler_method
    %v(i+1) = v(i) + g(x(i),v(i))*dt
    %x(i+1) = x(i) + f(x(i),v(i+1))*dt
    %standard Euler uses v(i) to compute x(i+1) instead of v(i+1)
    
    %given dt, forwardDyn.m assumes symplectic Euler
    %integrate qacc to update qvel
    qvel = qvel + qacc*dt;
end

%convert qvel to statedot
statedot = qvelToQdot(qvel, state(isorient), HT_world(1:3,1:3,1));

ydot = [statedot; qacc; u.err];

%additional outputs for data logging
out = OdeOutput();
out.HT_parent = HT_parent;
out.HT_world = HT_world;
out.contacts = contacts;
out.u = u.cmd;
out.vc = fdout.vc0;
% out.vc = fdout.vc; %DEBUGGING

out.fc = fdout.modelf_contact;
out.fact = fdout.modelf_act;
out.iter = fdout.iter;
out.cost = fdout.cost;



