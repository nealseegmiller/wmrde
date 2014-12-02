function [f,J] = odeWgc(p,vc,Rw,dz)
%Open Dynamics Engine wheel-ground contact model
%some modifications: fz dependent, rolling resistance

if nargout == 0
    %nothing to initialize
    return
end

scale = 1e3; %kN to N


Kp=p(1)*scale;
Kd=p(2)*scale;
%should be > 0, may be Inf
mux = p(3);
muy = p(4);
%should be <= 0
Cx = p(5)*scale; %= 1/force dependent slip
Cy = p(6)*scale;
Crr = p(7)*scale; %rolling resistance (NOT IN ODE!)


vx=vc(1,:);
vy=vc(2,:);
vz=vc(3,:);

nw=length(vx);

%compute normal force
fz = max(-Kp*dz - Kd*vz,0); %normal force must be positive
fz(dz > 0) = 0; %normal force must be zero if not in contact

fx = Cx*vx + Crr*Rw;
fy = Cy*vy;


%check limits
isover_x = abs(fx) > mux*fz;
isover_y = abs(fy) > muy*fz;

sign_fx_over = sign(fx(isover_x));
sign_fy_over = sign(fy(isover_y));

fx(isover_x) = sign_fx_over .* (mux*fz(isover_x));
fy(isover_y) = sign_fy_over .* (muy*fz(isover_y));

f=[fx;fy;fz];

if nargout > 1
    %Jacobian of [fx fy fz]' wrt [vx vy vz Rw dz]'
    J = zeros(3,5,nw);
    
    %d fz wrt dz,vz
    dfz_ddz = -Kp*ones(1,nw);
    dfz_dvz = -Kd*ones(1,nw);

    dfz_ddz(fz<=0) = 0;
    dfz_dvz(fz<=0) = 0;
    
    J(3,5,:) = dfz_ddz;
    J(3,3,:) = dfz_dvz;
   
    
    %d fx,fy wrt vx vy Rw 
    dfx_dvx = Cx*ones(1,nw);
    dfx_dRw = Crr*ones(1,nw);
    dfy_dvy = Cy*ones(1,nw);
    
    dfx_dvx(isover_x) = 0;
    dfx_dRw(isover_x) = 0;
    dfy_dvy(isover_y) = 0;

    J(1,1,:) = dfx_dvx;
    J(1,4,:) = dfx_dRw;
    J(2,2,:) = dfy_dvy;
    
    %d fx,fy wrt dz,vz
    dfx_dfz = zeros(1,nw);
    dfy_dfz = zeros(1,nw);
    
    dfx_dfz(isover_x) = sign_fx_over*mux;
    dfy_dfz(isover_y) = sign_fy_over*muy;
    
    J(1,3,:) = dfx_dfz.*dfz_dvz;
    J(1,5,:) = dfx_dfz.*dfz_ddz;
    
    J(2,3,:) = dfy_dfz.*dfz_dvz;
    J(2,5,:) = dfy_dfz.*dfz_ddz;

end

% blah=1;
