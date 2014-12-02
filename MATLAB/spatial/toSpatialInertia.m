function out = toSpatialInertia(m,c,I)

%from Featherstone's "mcI" function
%Neal's changes: 
%-removed planar
%-removed subfunctions

% mcI  rigid-body inertia <--> mass, CoM and rotational inertia.
% rbi=mcI(m,c,I) and [m,c,I]=mcI(rbi) convert between the spatial or planar
% inertia matrix of a rigid body (rbi) and its mass, centre of mass and
% rotational inertia about the centre of mass (m, c and I).  In the spatial
% case, c is 3x1, I is 3x3 and rbi is 6x6.  In the planar case, c is 2x1, I
% is a scalar and rbi is 3x3.  In both cases, m is a scalar.  When c is an
% argument, it can be either a row or a column vector.  If only one
% argument is supplied then it is assumed to be rbi; and if it is a 6x6
% matrix then it is assumed to be spatial.  Otherwise, three arguments must
% be supplied, and if length(c)==3 then mcI calculates a spatial inertia
% matrix.  NOTE: (1) mcI(rbi) requires rbi to have nonzero mass; (2) if |c|
% is much larger than the radius of gyration, or the dimensions of the
% inertia ellipsoid, then extracting I from rbi is numerically
% ill-conditioned.

C = skew3(c);
out = [ I + m*(C*C'), m*C; m*C', m*eye(3) ];