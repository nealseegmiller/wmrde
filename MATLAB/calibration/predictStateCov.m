function P = predictStateCov(P0,Q,orient,R_body_to_world,vel,dt)
%INPUTS
%P0:    (ns x ns) state covariance at time i
%Q:     (nv x nv) spectral density of qvel (/dt?)
%orient:            orientation vector
%R_body_to_world:   3x3 Rotation from body to world coords, redundant with orientation! assumed to be consistent
%vel:   6x1 spatial velocity of body frame in body coords [wx wy wz vx vy vz]' 
%dt:    time i+1 - time i
%OUTPUT
%P at time i+1


ns = size(P0,1);
nv = size(Q,1);
olen = length(orient);

nf = ns-olen-3+1; %number of frames in WmrModel
[isorient,ispos,isjd] = stateIs(ns,nf);

[~,T_vel_to_orate,dTdo] = velToOrientrate(orient);
[~,T_orate_to_vel] = orientrateToVel(orient);
    
G = zeros(ns,nv);


temp1 = -R_body_to_world*skew3(vel(4:6))*T_orate_to_vel;
temp2 = zeros(olen);
for i=1:olen
    temp2(:,i) = dTdo(:,:,i)*vel(1:3); 
end

if 1    %discrete time
    G(isorient,1:3) = T_vel_to_orate*dt;
    G(ispos,4:6) = R_body_to_world*dt;
    G(isjd,7:end) = eye(nf-1)*dt;
    
    F = eye(ns);
    F(ispos,isorient) = temp1*dt;
    F(isorient,isorient) = F(isorient,isorient) + temp2*dt;

    P = F*P0*F' + G*Q*G';
else    %continuous time
    G(isorient,1:3) = T_vel_to_orate;
    G(ispos,4:6) = R_body_to_world;
    G(isjd,7:end) = eye(nf-1);
    
    F = zeros(ns);
    F(ispos,isorient) = temp1;
    F(isorient,isorient) = F(isorient,isorient) + temp2;

%     P = P0 + (F*P0 + P0*F' + G*(Q*dt)*G')*dt;
    P = P0 + (F*P0 + P0*F' + F*P0*F'*dt + G*(Q*dt)*G')*dt; %Stengel 1994, equation 4.2-75
end


