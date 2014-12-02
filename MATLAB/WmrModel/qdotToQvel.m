function qvel = qdotToQvel(qdot,orient,R_body_to_world)
%inverse of qdotToQvel()

if nargin < 3
    R_body_to_world = orientToRot(orient);
end

[ns,m]=size(qdot);
olen = length(orient);
% olen = SIZEORIENT(); %DEBUGGING
nf = ns - (olen + 3) + 1;
nv = nf + 5; %size of joint space vel


[isorient,ispos,isjd] = stateIs(ns);

qvel=zeros(nv,m);
qvel(7:end,:) = qdot(isjd,:);
qvel(4:6,:) = R_body_to_world'*qdot(ispos,:);
qvel(1:3,:) = orientrateToVel(orient,qdot(isorient,:));


