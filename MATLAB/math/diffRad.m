function out = diffRad(theta1,theta2)
%difference angles in radians, output range of -pi to pi
%out = theta1-theta2
%theta1, theta2 are arrays of angles of the same size

max = pi;
out = wrapRad(theta1-theta2,max);
