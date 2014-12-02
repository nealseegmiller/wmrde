function vslip = universalBsm(c,v,grav,Is,R_body_to_world,lon,lat,ang)

%universal body slip model
%INPUTS
%c:     nx1 coefficients for slip parameterization
%v:     6x1 spatial velocity of body frame (in body coords) least-squares slip solution
%grav:  scalar, magnitude of gravitational acceleration, e.g. 9.81 m/s^2
%Is:    6x6 spatial inertia of WMR (in body coords)
%R_body_to_world:  3x3 rotation matrix, body to world coords
%(optional)
%lon,lat,ang 3x1 unit vectors in slip directions (in body coords). Default aligned with body frame x,y,z axes
%OUTPUTS
%vslip, 6x1, spatial slip velocity of body frame (in body coords)

%lon,lat,ang directions
if nargin < 6
    %aligned with body frame axes
    lon = [1 0 0]';
    lat = [0 1 0]';
    ang = [0 0 1]';
end

%for indexing spatial vectors
ia = 1:3; %angular part
il = 4:6; %linear part

%spatial gravity acceleration (in body coords)
g = zeros(6,1);
% g(il) = HT_body_to_world(1:3,1:3)'*[0 0 -grav];
g(il) = -R_body_to_world(3,1:3)'*grav; %faster!

%compute inertial force, g-a
f = Is*g - crossForce(v)*Is*v; %(signs verified)


%transform to cg frame
%Plucker transforms from body to cg
[~,rcg,~] = fromSpatialInertia(Is);

% HT_cg_to_body = eye(4); HT_cg_to_body(1:3,4) = rcg;
% Xm = invHTToPlucker(HT_cg_to_body); 
% Xf = HTToPlucker(HT_cg_to_body)';

%faster
tmp = -skew3(rcg);
Xm = [eye(3), zeros(3); tmp, eye(3)];
Xf = [eye(3), tmp; zeros(3), eye(3)];
    
vcg = Xm*v;
fcg = Xf*f; %moment goes to zero

%DEBUGGING, scale up flon,flat so all c coefficients are about the same size
scale=10;

%precompute
flon=fcg(il)'*lon*scale;
flat=fcg(il)'*lat*scale;
fz = fcg(il)'*ang; %normal force

vlon=vcg(il)'*lon;
vang=vcg(ia)'*ang;

vcg_slip = zeros(6,1); %spatial slip velocity of cg

vcg_slip(il) = (c(1)*flon/fz*vlon + c(2)*vlon)*lon + (c(3)*flat/fz*vlon)*lat; %lon/lat slip
vcg_slip(ia) = (c(4)*flat/fz*vlon + c(5)*vlon + c(6)*vang)*ang; %ang slip


%transform back to body frame
vslip = Xf'*vcg_slip; %Xf' = inv(Xm)






