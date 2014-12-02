function M=quatmatr(q)
%see quatmatl

w=q(1); x=q(2); y=q(3); z=q(4);

M = [w -x -y -z; 
     x  w  z -y; 
     y -z  w  x; 
     z  y -x  w];
 
