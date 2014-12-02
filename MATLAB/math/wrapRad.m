function theta = wrapRad(theta,mymax)
%wrap angles in radians to range of 2*pi
%INPUTS
%theta, array of angles
%mymax, max angle value (default pi)

if nargin < 2
    mymax=pi;
end
myrange = 2*pi;

theta = mod(theta + mymax, myrange) + (mymax - myrange); 