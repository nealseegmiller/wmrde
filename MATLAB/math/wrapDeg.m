function theta = wrapDeg(theta,mymax)
%wrap angles in degrees to range of 360
%INPUTS
%theta, array of angles
%mymax, max angle value (default 180)

if nargin < 2
    mymax=180;
end
myrange = 360;

theta = mod(theta + mymax, myrange) + (mymax - myrange); 