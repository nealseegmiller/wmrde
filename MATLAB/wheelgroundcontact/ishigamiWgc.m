function [f,J] = ishigamiWgc(p,vc,Rw,dz)
%ref. Ishigami et al., Terramechanics-based model for steering maneuver of planetary exploration rovers on loose soil
%also ref. Jia et al., Fast computation of wheel-soil interactions ... ICRA 2011
%NOTE
%-non-zero fx for s==0, rolling-resistance


if nargin==1
    %nothing to initialize
    return
end

vx=vc(1,:);
vy=vc(2,:);
vz=vc(3,:);
nw = length(vx);


f = subfn(p,vc,Rw,dz);

if nargout > 1
    J = zeros(3,5,nw);

    myeps = 1e-6;
    for i = [1 2 4 5] %finite diff, all except vz
        if i==1
            fp = subfn(p,[vx+myeps;vy;vz],Rw,dz);
        elseif i==2
            fp = subfn(p,[vx;vy+myeps;vz],Rw,dz);
        elseif i==4
            fp = subfn(p,vc,Rw+myeps,dz);
        elseif i==5
            fp = subfn(p,vc,Rw,dz+myeps);
        end      
        
        %forward difference
        J(:,i,:) = reshape((fp-f)/myeps,[3 1 nw]);

    end
    
    Kd = p(1);
    J(3,3,:) = -Kd; %vz

    fz=f(3,:);
    J(3,:,fz==0)=0;
    
end

end

function f = subfn(p,vc,Rw,dz)

vx=vc(1,:);
vy=vc(2,:);
vz=vc(3,:);

Kd = p(1); %(N/(m/s)) damping coefficient, fz -= Kd*vz, necessary?
r = p(2); %(m) wheel radius
b = p(3); %(m) wheel width
c = p(4); %(kPa) cohesion stress
phi = p(5); %(rad) friction angle
Xc = p(6); %(rad) soil distractive angle, ~pi/4-phi/2;
kc = p(7); %(N/m^2) pressure sinkage module
kphi = p(8); %(N/m^3) pressure sinkage module
n = p(9); %sinkage exponent
%a0,a1 determine theta m (max normal pressure)
a0 = p(10); 
a1 = p(11);
rhod = p(12); %(N/m^3) soil density, not kg/m^3 in paper must be a typo. units don't work out
lambda = p(13); %wheel sinkage ratio
if length(p) > 14
    %soil deformation module
    kx_m = p(14); %(m)
    kx_b = p(15); %(m)
    ky_m = p(16); %(m)
    ky_b = p(17); %(m)
    isotropic=false;
else
    K = p(14); %(m) soil deformation module
    isotropic=true;
end


np = 61; %number of points, for discretization of stress distributions
nw = length(vx);

%convert to slip ratio and angle
[s,beta] = calcSlip(vx,vy,Rw,2);


%entry and exit angles
h = max(-dz,0);
theta_f = acos(1-h/r);
theta_r = -acos(1-(lambda/r)*h);

%discretize angles
theta = zeros(np,nw);
% for i = 1:nw
for i = find(dz<0)
    theta(:,i) = theta_r(i):(theta_f(i)-theta_r(i))/(np-1):theta_f(i);
end
%precompute:
ctheta = cos(theta);
stheta = sin(theta);

%normal stress distribution
sigma = zeros(np,nw);
theta_m = (a0 + a1*s).*theta_f; %okay, intuitively wheel could rise out of ground when braking

temp1 = r^n*(kc/b + kphi);
temp2 = (theta_f - theta_m)./(theta_m - theta_r); %1xnw
for i = 1:nw
    isfront = theta(:,i) >= theta_m(i);
    sigma(isfront,i) = temp1*(ctheta(isfront,i) - cos(theta_f(i))).^n;
    sigma(~isfront,i) = temp1*(cos(theta_f(i) - (theta(~isfront,i) - theta_r(i))*temp2(i)) - cos(theta_f(i))).^n;
end

%soil deformations
jx = zeros(np,nw);
jy = zeros(np,nw);

taux = zeros(np,nw); 
tauy = zeros(np,nw);

for i = 1:nw

    
    %these are based on s = (Rw-Vx)/Rw
    %vx(theta) = Rw - Vx*cos(theta)
    %   sub Vx = Rw*(1-s)
    %vx(theta) = Rw - Rw*(1-s)*cos(theta) = Rw*(1-(1-s)*cos(theta))
    %vy(theta) = -Vx*tan(beta)
    %vy(theta) = -Rw*(1-s)*tan(beta)
    %large deformation okay because of exponential function in shear eqn
    
    jx(:,i) = r*( theta_f(i) - theta(:,i) - (1-s(i))*(sin(theta_f(i)) - stheta(:,i)) );
%     jy(:,i) = -r*(1-s(i))*(theta_f(i)-theta(:,i))*tan(beta(i)); %negate!
    
    jy(:,i) = -r*(1-abs(s(i)))*(theta_f(i)-theta(:,i))*tan(beta(i)); %DEBUGGING, make symmetric

    if isotropic
        %Jia-isotropic, direction via shear displacement vs. velocity
        jnorm = sqrt(jx(:,i).^2 + jy(:,i).^2);
        tau = (c+sigma(:,i)*tan(phi)).*(1-exp(-jnorm/K));

        taux(:,i) = tau.*jx(:,i)./jnorm;
        tauy(:,i) = tau.*jy(:,i)./jnorm;
    else
        %Ishigami-anisotropic

        kx=kx_m*abs(beta(i))+kx_b;
        ky=ky_m*abs(beta(i))+ky_b;

        taux(:,i) = (c+sigma(:,i)*tan(phi)).*sign(jx(:,i)).*(1-exp(-abs(jx(:,i))/kx));
        tauy(:,i) = (c+sigma(:,i)*tan(phi)).*sign(jy(:,i)).*(1-exp(-abs(jy(:,i))/ky));
        
    end

end

taux(isnan(taux))=0;
tauy(isnan(tauy))=0;


if 0
    %compute bulldozing reaction force as f(theta)   
    htheta = r*(ctheta - ones(np,1)*cos(theta_f)); %sinkage h as a function of theta

    D1 = cot(Xc) + tan(Xc+phi);
    D2 = cot(Xc) + cot(Xc)^2/cot(phi);
    Rb = D1*(c*htheta + (D2*rhod/2)*htheta.^2);
    Fs = Rb.*(r-htheta.*ctheta);
    
    %compute bulldozing force coefficient for each wheel
    %no dependence of magnitude on beta is reasonable if quasi-static assumption
    %but beta must determine direction
    bfc = -sign(beta); %discontinuity in derivative wrt beta at beta=0
%     bfc = -sign(beta).*(1-exp(-abs(beta)/2e-2));
    
    Fs = (ones(np,1)*bfc) .* Fs;
else
    Fs = zeros(np,nw); %DEBUGGING
end

%integrate to compute forces
fx = zeros(1,nw);
fy = zeros(1,nw);
fz = zeros(1,nw);
for i = 1:nw
    fx(i) = r*b*trapz(theta(:,i), taux(:,i).*ctheta(:,i) - sigma(:,i).*stheta(:,i));
    fy(i) = trapz(theta(:,i), r*b*tauy(:,i) + Fs(:,i));
    fz(i) = r*b*trapz(theta(:,i), taux(:,i).*stheta(:,i) + sigma(:,i).*ctheta(:,i));
end

%add damping to normal force
fz = fz - Kd*vz;

%normal force must be positive
fz(fz<0)=0;

if 1
    %reverse sign for reverse, TODO check this!
    Vx=vx+Rw;
    isrev = Rw<0;
    fx(isrev) = -fx(isrev);
    isrev = Vx<0;
    fy(isrev) = -fy(isrev);
end


f=[fx;fy;fz];

end









