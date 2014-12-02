function M=quatmatl(q)
%convert to matrix for quaternion multiplication
%http://en.wikipedia.org/wiki/Rotation_formalisms_in_three_dimensions#Quaternions (different order convention)
%http://www.cprogramming.com/tutorial/3d/quaternions.html

% q3 = q1 (x) q2
% q3 = quatmatl(q1)*q2;
% q3 = quatmatr(q2)*q1;

w=q(1); x=q(2); y=q(3); z=q(4);

M = [w -x -y -z; 
     x  w -z  y; 
     y  z  w -x; 
     z -y  x  w];

 