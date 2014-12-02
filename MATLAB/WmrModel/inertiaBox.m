function I = inertiaBox(m,x,y,z)
%calculate the 3x3 rotational inertia of a rectangular box about it's center of mass 
%given: mass, length in the x,y,z dimensions

Ixx = m*(y^2 + z^2)/12;
Iyy = m*(x^2 + z^2)/12;
Izz = m*(x^2 + y^2)/12;

I = diag([Ixx, Iyy, Izz]);