function  vcross = crossMotion( v )

%from Featherstone's "crm" function

%Neal's changes:
%-removed planar code

% crm  spatial/planar cross-product operator (motion).
% crm(v)  calculates the 6x6 (or 3x3) matrix such that the expression
% crm(v)*m is the cross product of the motion vectors v and m.  If
% length(v)==6 then it is taken to be a spatial vector, and the return
% value is a 6x6 matrix. 


vcross = [0    -v(3)  v(2)   0     0     0   ;
	      v(3)  0    -v(1)   0     0     0   ;
	     -v(2)  v(1)  0      0     0     0   ;
	      0    -v(6)  v(5)   0    -v(3)  v(2);
	      v(6)  0    -v(4)   v(3)  0    -v(1);
	     -v(5)  v(4)  0     -v(2)  v(1)  0 ];

