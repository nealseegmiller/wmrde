function [rol,pit] = inclinationToEuler(inc_x,inc_y)
%convert inclination angles to Euler roll and pitch angles
%the difference is negligible for small angles

%Convention: R^(world)_(body) = Rotz(yaw)*Roty(pitch)*Rotx(roll)
%   body frame z axis in world coords = [sp*cr, -sr, cp*cr]'
%   inc_x = atan(-z(2)/z(3)), inclination angle about x axis
%   inc_y = atan( z(1)/z(3)), inclination angle about y axis

pit = inc_y;
rol = atan(cos(pit).*tan(inc_x));