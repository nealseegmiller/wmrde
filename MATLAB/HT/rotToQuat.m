function q = rotToQuat(R)
%convert rotation matrix to quaternion
%http://en.wikipedia.org/wiki/Rotation_matrix#Quaternion

Qxx=R(1,1);
Qyy=R(2,2);
Qzz=R(3,3);
Qzy=R(3,2); Qyz=R(2,3);
Qxz=R(1,3); Qzx=R(3,1);
Qyx=R(2,1); Qxy=R(1,2);

t = Qxx+Qyy+Qzz;
r = sqrt(1+t);
w = 0.5*r;

if t > 0
    s = 0.5/r;
    x = (Qzy-Qyz)*s;
    y = (Qxz-Qzx)*s;
    z = (Qyx-Qxy)*s;
else
    %avoid singularity
    x = sign(Qzy-Qyz)*0.5*sqrt(1+Qxx-Qyy-Qzz);
    y = sign(Qxz-Qzx)*0.5*sqrt(1-Qxx+Qyy-Qzz);
    z = sign(Qyx-Qxy)*0.5*sqrt(1-Qxx-Qyy+Qzz);
end

q = [w x y z]';
