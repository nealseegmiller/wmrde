function R = quatToRot(q)
%convert quaternion to rotation matrix
%http://en.wikipedia.org/wiki/Rotation_matrix#Quaternion
%valid for non-unit quaternions

q = q/sqrt(sum(q.^2)); %normalize
w=q(1); x=q(2); y=q(3); z=q(4);

wx = w*x; wy = w*y; wz = w*z;
xx = x*x; xy = x*y; xz = x*z;
yy = y*y; yz = y*z; 
zz = z*z;

R=[...
1-2*(yy+zz),   2*(xy-wz),   2*(xz+wy)
  2*(xy+wz), 1-2*(xx+zz),   2*(yz-wx)
  2*(xz-wy),   2*(yz+wx), 1-2*(xx+yy)
];



