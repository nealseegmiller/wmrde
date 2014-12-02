function [eulerrate,T,dTdeuler] = velToEulerrate(euler,angvel)
%convert angular velocity (in body coords) to Euler angle rates
%singularity at pit = pi/2
%INPUTS
%euler:     3 x 1 euler angles
%angvel:    (optional) 3 x 1, angular velocity
%OUTPUTS
%eulerrate: 3 x 1, d/dt euler. [] if no angvel input
%T:         3 x 3, transform eulerrate = T*angvel
%dTdeuler:  3 x 3 x 3, dT/deuler


rol=euler(1);
pit=euler(2);
%yaw not required!

sR = sin(rol);
cR = cos(rol);
cP = cos(pit);
tP = tan(pit);

%ASSUME 321 ROTATION ORDER! SAE CONVENTION
T = [1, sR*tP, cR*tP
     0,    cR,   -sR
     0, sR/cP, cR/cP];

eulerrate = [];
if nargin > 1
    eulerrate = T*angvel;
end

if nargout > 2

    %derivative of T wrt Euler angles
    dTdrol = [0, cR*tP, -sR*tP
              0,   -sR,    -cR
              0, cR/cP, -sR/cP];

    sP=sin(pit);
    cP2=cP^2;
    dTdpit = [0,    sR/cP2,    cR/cP2
              0,         0,         0
              0, sR*sP/cP2, cR*sP/cP2];

    dTdeuler = cat(3, dTdrol, dTdpit, zeros(3));

end

