function [f,J] = pacejka_wgc(p,vc,Rw,dz)
%Pacejka (or BNP) wheel-ground contact model
%ref. VCRware model in Brach & Brach 09, Tire Models for Vehicle Dynamic Simulation and Accident Reconstruction

if nargin==1
    %nothing to initialize
    return
end

Kp = p(1);
Kd = p(2);

mu = p(3);

Bs = p(4);
Cs = p(5);
Ds = p(6);
Es = p(7);
Ks = p(8);

Ba = p(9);
Ca = p(10);
Da = p(11);
Ea = p(12);
Ka = p(13);

Crr = p(14); %rolling resistance

vx=vc(1,:);
vy=vc(2,:);
vz=vc(3,:);

nw = length(vx); %number of wheels

%compute normal force
fz = max(-Kp*dz - Kd*vz,0); %normal force must be positive
fz(dz > 0) = 0; %normal force must be zero if not in contact

%convert to slip ratio and angle
[s,alpha,ds_dvx,ds_dRw,dalpha_dvx,dalpha_dvy,dalpha_dRw] = calcSlip(vx,vy,Rw,1);


%Nicolas Comstock singularity at s==0 or alpha==0
%TODO, a better way? must choose eps carefully
%affects Jacobian accuracy
eps = 1e-10;
s(s >= 0 & s < eps) = eps;
s(s < 0 & s > -eps) = -eps;
alpha(alpha >= 0 & alpha < eps) = eps;
alpha(alpha < 0 & alpha > -eps) = -eps;


%Pacejka, forces for only lon or only lat slip
Fxs = Ds*sin(Cs*atan(Bs*(1-Es)*Ks*s + Es*atan(Bs*Ks*s)));
temp1 = (2/pi)*alpha;
Fya = Da*sin(Ca*atan(Ba*(1-Ea)*Ka*temp1 + Ea*atan(Ba*Ka*temp1)));

%handle reverse
Vx = vx+Rw;
sgn = ones(1,nw); 
sgn(Vx<0)=-1;
    
if 1
    %Nicolas Comstock, forces combined lon/lat slip
    %See eqn 9 and 10 in Brach 2009
    s2 = s.*s;
    Fxs2 = Fxs.*Fxs;
    Fya2 = Fya.*Fya;
    tan_alpha = tan(alpha);
    tan_alpha2 = tan_alpha.*tan_alpha;
    temp2 = abs(Fxs.*Fya) ./ sqrt( s2.*Fya2 + Fxs2.*tan_alpha2 );

    dfx_dfz = sgn .* temp2 .* s * mu;
    dfy_dfz = -sgn .* temp2 .* tan_alpha * mu; %I added the negative sign, Brach defines different coordinate system for force (Fig 2)

else
    %DEBUGGING, Pacejka only
    dfx_dfz = sgn.*Fxs*mu;
    dfy_dfz = -sgn.*Fya*mu;
end
dfx_dfz = dfx_dfz + Crr*Rw; %add rolling resistance

fx = dfx_dfz .* fz;
fy = dfy_dfz .* fz;

f=[fx;fy;fz];

if nargout > 1
    
    J = zeros(3,5,nw); %J(:,:,i) is the Jacobian of [fx fy fz]' wrt [vx vy vz Rw dz]' for wheel i
    
    %d fz wrt dz,vz
    dfz_ddz = -Kp*ones(1,nw);
    dfz_dvz = -Kd*ones(1,nw);
    dfz_ddz(fz<=0) = 0;
    dfz_dvz(fz<=0) = 0;
    
    J(3,5,:) = dfz_ddz;
    J(3,3,:) = dfz_dvz;
    
    %derivatives of fx,fy wrt abs(s),abs(alpha)
    [dfx_ds,dfy_ds,dfx_dalpha,dfy_dalpha] = pacejkaDerivatives(Ba,Bs,Ca,Cs,Da,Ds,Ea,Es,Ka,Ks,alpha,fz,mu,s,sgn);
    
    
    J(1,1,:) = dfx_ds.*ds_dvx + dfx_dalpha.*dalpha_dvx;
    J(1,2,:) = dfx_dalpha.*dalpha_dvy;
    J(1,3,:) = dfx_dfz.*dfz_dvz;
%     J(1,4,:) = dfx_ds.*ds_dRw + dfx_dalpha.*dalpha_dRw;
    J(1,4,:) = dfx_ds.*ds_dRw + dfx_dalpha.*dalpha_dRw + Crr*fz; %add rolling resistance
    J(1,5,:) = dfx_dfz.*dfz_ddz;
    
    J(2,1,:) = dfy_ds.*ds_dvx + dfy_dalpha.*dalpha_dvx;
    J(2,2,:) = dfy_dalpha.*dalpha_dvy;
    J(2,3,:) = dfy_dfz.*dfz_dvz;
    J(2,4,:) = dfy_ds.*ds_dRw + dfy_dalpha.*dalpha_dRw;
    J(2,5,:) = dfy_dfz.*dfz_ddz;
    
end

% return


