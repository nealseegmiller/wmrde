function Q = bodyVelNoise(p,v)
%INPUT
%p:     (4x1) parameter vector [lon lat ang att]
%v:     (6x1) spatial velocity of body frame in body coords [wx wy wz vx vy vz]' 
%OUTPUT
%Q:     (6x6 matrix) spectral density of body frame velocity noise

%TODO, make alternative function that converts wheel-level to body-level noise using wheel Jacobian

% if norm(v(4:6)) > 1e-3 %only add noise if moving

    Q_std = zeros(6,1);
    Q_std([4 5 3 1 2]) = p([1 2 3 4 4]); %[vx vy wz wx wy] = [lon lat ang att att]
    Q = diag(Q_std.^2);

% end



