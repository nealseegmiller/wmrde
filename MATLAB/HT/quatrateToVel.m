function [angvel,T] = quatrateToVel(q,quatrate)
%convert quaternion rate to angular velocity (in body coords)

T = quatmatl(invertQuat(q));
T = 2*T(2:4,:);

angvel = [];
if nargin > 1
    angvel=T*quatrate;
end
