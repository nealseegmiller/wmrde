function p = pointsOnCircle(rad,nv,close,ax)
%INPUTS
%rad:   circle radius
%nv:    number of vertices
%close: if true, make a closed polygon with nv+1 points (the last point is the same as the first), else nv points
%ax:    (optional) axis of symmetry 1,2,3
%OUTPUT
%p:     a list of points equally distributed along the circumference of a circle, clockwise ordered, origin at (0,0)
%       2xnp if no ax input
%       else 3xnp with p(ax,:) = 0

if close
    np = nv+1;
    angles = linspace(0,2*pi,np);
else
    np = nv;
    angles = linspace(0,2*pi-2*pi/np,np);
end
p_(1,:) = rad*cos(angles);
p_(2,:) = -rad*sin(angles);
%clockwise ordered vertices, required for polybool

if nargin < 4
    p = p_;
else
    p = zeros(3,np);
    p(1:3~=ax,:) = p_;   
end

