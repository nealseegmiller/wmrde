function [m, c, I] = fromSpatialInertia(rbi)
%convert from spatial inertia to mass, center of mass, and moment of inertia
%from Featherstone's "mcI" function

m = rbi(6,6);
mC = rbi(1:3,4:6);
c = unskew3(mC)/m;
I = rbi(1:3,1:3) - mC*mC'/m;