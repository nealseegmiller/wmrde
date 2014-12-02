function cp = contactAnglesToPoints(radius, angles)
%compute contact points in wheel coords given contact angles
%convention:
%contact angle=0 corresponds to cp=[0,0,-rad]'
%ccw about y axis is positive
%INPUTS
%radius:    scalar
%angles:	1 x np

np = length(angles);
cp = zeros(3,np);
cp(1,:) = -sin(angles)*radius;
cp(3,:) = -cos(angles)*radius;


