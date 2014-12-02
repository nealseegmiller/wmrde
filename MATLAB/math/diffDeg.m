function out = diffDeg(theta1,theta2)
%difference angles in degrees, output range of -180 to 180
%out = theta1-theta2
%theta1, theta2 are arrays of angles of the same size

max = 180;
out = wrapDeg(theta1-theta2,max);
