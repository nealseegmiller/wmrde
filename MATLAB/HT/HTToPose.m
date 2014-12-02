function [orient,pos] = HTToPose(HT,usequaternion)
%get pose vector from matrix

R = HT(1:3,1:3);
pos = HT(1:3,4);

if usequaternion
    orient = rotToQuat(R);
else
    %use Euler angles
    orient = rotToEuler(R);
end

