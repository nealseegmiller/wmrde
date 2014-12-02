function out = HTToPlucker(in)

%converts a 4x4 homogeneous transform to 6x6 Plucker transform
%from Featherstone's "pluho" function

%the following identity is used:
%-E*skew(r) == skew(-E*r)*E
%-E*skew(r)*E' == skew(-E*r)
% E*skew(r)*E' == skew(E*r)


% pluho  convert Plucker <--> 4x4 homogeneous coordinate transform.
% X=pluho(T) and T=pluho(X) convert between a Plucker coordinate transform
% matrix X and a 4x4 homogeneous coordinate transform matrix T.  If the
% argument is a 6x6 matrix then it is taken to be X, otherwise T.
% NOTE: the 4x4 matrices used in 3D graphics (e.g. OpenGL and Matlab handle
% graphics) are displacement operators, which are the inverses of
% coordinate transforms.  For example, to set the 'matrix' property of a
% Matlab hgtransform graphics object so as to rotate its children by an
% angle theta about the X axis, use inv(pluho(rotx(theta))).

% Formulae:
%   X = [  E   0 ]   T = [ E  -Er ]
%       [ -Erx E ]       [ 0   1  ]

E = in(1:3,1:3);
mEr = in(1:3,4);			% - E r
out = [ E, zeros(3); skew3(mEr)*E, E ];