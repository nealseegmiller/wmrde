function out = invHTToPlucker(in)
%converts the inverse of a 4x4 homogeneous transform to a 6x6 Plucker transform
%faster than HTToPlucker(invertHT(in)), but still slow!

%ht = [E, -Er   inv(ht) = [E',-E'*(-Er)
%      0,  1]              0,  1]

% E = in(1:3,1:3);
% mEr = in(1:3,4);      %-E r
% 
% %invert
% E_ = E';
% mEr_ = -E'*mEr;
% 
% out = [ E_, zeros(3); skew3(mEr_)*E_, E_ ];


%this is faster
%   E_ = E';
%   skew3(mEr_)*E_ = skew3(-E'*mEr)*E'
%                  =-skew3(E'*mEr)*E'
%                  =-(E'*skew3(mEr)*E)*E' 
%                  =-E'*skew3(mEr)


Et = in(1:3,1:3)'; %E'
mEr = in(1:3,4);
% out = [Et, zeros(3); -Et*skew3(mEr), Et];

%faster?
out = zeros(6);
out(1:3,1:3) = Et;
out(4:6,4:6) = Et;
out(4:6,1:3) = -Et*skew3(mEr);



