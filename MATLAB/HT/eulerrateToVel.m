function [angvel,T] = eulerrateToVel(euler,eulerrate)
%convert Euler angle rates to angular velocity (in body coords)

rol = euler(1);
pit = euler(2);
%yaw not required!

sR = sin(rol); cR = cos(rol);
sP = sin(pit); cP = cos(pit);

%ASSUME 321 ROTATION ORDER! THIS IS THE SAE CONVENTION
T = [1,   0,    -sP
     0,  cR,  sR*cP
     0, -sR,  cR*cP];

angvel = [];
if nargin > 1
    angvel = T*eulerrate; %in body coords!
end

