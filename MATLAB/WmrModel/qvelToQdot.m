function qdot = qvelToQdot(qvel,orient,R_body_to_world)

%Explanation of qvel vs. qdot
%the first elements of state (q) are pose of the body frame in world coords,
%orientation:
%[rol pit yaw] (if using Euler angles)
%[q1 q2 q3 q4] (if using quaternion)
%position: [x y z]
%the remaining elements are joint displacements

%qvel differs from d/dt(state) (qdot) in that the derivative of pose wrt time is replaced with [wx wy wz vx vy vz]' (in body coords)
%The first 3 elements are angular velocities, the next 3 are linear velocities 
%Conversion to qdot requires transforming angular velocity to Euler angle or quaternion rates

%Likewise, qacc differs from d^2/dt^2 state (qddot) in that the 2nd derivative of pose wrt time is replaced with spatial acceleration

%INPUTS
%qvel:              nv x 1 vector (or nv x m matrix)
%orient:            3 x 1 vector of Euler angles or 4 x 1 quaternion
%R_body_to_world:   3 x 3 rotation matrix, redundant with orientation (assumed to be consistent)

if nargin < 3
    R_body_to_world = orientToRot(orient);
end

[nv,m] = size(qvel);
olen = length(orient);
% olen = SIZEORIENT(); %DEBUGGING
nf = nv - 5;
ns = nf - 1 + olen + 3; %number of states

[isorient,ispos,isjd] = stateIs(ns);

qdot = zeros(ns,m);
qdot(isjd,:) = qvel(7:end,:);
qdot(ispos,:) = R_body_to_world*qvel(4:6,:);
qdot(isorient,:) = velToOrientrate(orient,qvel(1:3,:));


