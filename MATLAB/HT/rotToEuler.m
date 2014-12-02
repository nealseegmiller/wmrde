function euler = rotToEuler(R)

%ASSUME 321 ROTATION ORDER! THIS IS THE SAE CONVENTION
%atan2 returns values in the range [-pi,pi]
yaw = atan2(R(2,1), R(1,1)); 
rol = atan2(R(3,2), R(3,3));
srol=sin(rol); crol=cos(rol);
if abs(crol) > abs(srol)
    cpit = R(3,3)/crol; %singularity at rol = +/- pi/2
else
    cpit = R(3,2)/srol;
end
pit = atan2(-R(3,1), cpit);
euler = [rol; pit; yaw];
