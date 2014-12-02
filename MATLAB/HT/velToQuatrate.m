function [qrate,T,dTdq] = velToQuatrate(q,angvel)
%convert angular velocity (in body coords) to quaternion rate
%http://en.wikipedia.org/wiki/Rotation_formalisms_in_three_dimensions#Quaternion_.E2.86.94_angular_velocities
%INPUTS
%q:         4 x 1 quaternion
%angvel:    (optional) 3 x 1, angular velocity
%OUTPUTS
%quatrate:  4 x 1, d/dt q. [] if no angvel input
%T:         4 x 3, transform quatrate = T*angvel
%dTdq:      4 x 3 x 4, dT/dq

T = .5*quatmatl(q);
T = T(:,2:4);

qrate = [];
if nargin > 1
    qrate=T*angvel;
end

if nargout > 2
    
    %derivative of T wrt quaternion
    %computed using symbolic toolbox
    dTdq = zeros(4,3,4);
    
%     T =
%     [ -q2/2, -q3/2, -q4/2]
%     [  q1/2, -q4/2,  q3/2]
%     [  q4/2,  q1/2, -q2/2]
%     [ -q3/2,  q2/2,  q1/2]

    dTdq(:,:,1) = [ 0  0  0; 
                    1  0  0; 
                    0  1  0; 
                    0  0  1]*.5;

    dTdq(:,:,2) = [-1  0  0; 
                    0  0  0; 
                    0  0 -1; 
                    0  1  0]*.5;

    dTdq(:,:,3) = [ 0 -1  0; 
                    0  0  1; 
                    0  0  0; 
                   -1  0  0]*.5;

    dTdq(:,:,4) = [ 0  0 -1; 
                    0 -1  0; 
                    1  0  0; 
                    0  0  0]*.5;
end
