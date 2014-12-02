function I = inertiaCylinder(m,r,h,ax)
%calculate the 3x3 rotational inertia of a cylinder about its center of mass
%given: mass, radius, height, and axis of symmetry

Idiag = m*(3*r^2 + h^2)/12*ones(3,1);
Idiag(ax) = m*r^2/2;

I = diag(Idiag);


