function [ydot,out] = odeKin(time,y,mdl,surfaces,contacts,pathno)
%ordinary differential equation for kinematic simulation
%INPUTS
%time:      scalar
%y:         state
%mdl:       WmrModel object
%surfaces:  cell array of surface objects for terrain
%contacts:  1 x nw array of WheelContactGeom objects, or
%           1 x nt array of TrackContactGeom objects
%OUTPUTS
%ydot:      d/dt y
%out:       OdeOutput object

if nargin < 6
    pathno = [];
end

%get from WmrModel object
nf = mdl.nf;
nv = nf+5;

isorient = stateIs(length(y));

state=y;

%convert state to homogeneous transforms
[HT_parent, HT_world] = stateToHT(mdl,state);

%update contact geometry
contacts = updateModelContactGeom(mdl, surfaces, HT_world, mdl.min_npic, contacts);

%get control inputs
if isempty(pathno)
    u_ = feval(mdl.controller_fh,mdl,time,state);
else
    u_ = feval(mdl.controller_fh,mdl,time,state,pathno);
end
u = NaN(nv,1);
u(mdl.actframeinds+5) = u_;

tstart = tic;

%compute velocity, DAE
[qvel,vc] = forwardVelKin(mdl,state,u,HT_world,contacts);
    
motiontime = toc(tstart);

statedot = qvelToQdot(qvel, state(isorient), HT_world(1:3,1:3,1)); %convert qvel to statedot

ydot=statedot;

%additional outputs for data logging
out = OdeOutput();
out.HT_parent = HT_parent;
out.HT_world = HT_world;
out.contacts = contacts;
out.u = u_;
out.vc = vc;
out.motiontime = motiontime;



