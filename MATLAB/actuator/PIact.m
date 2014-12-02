function [tau,err,dtau_du] = PIact(p,ucmd,u,interr)
%actuator model with proportional integral (PI) velocity control

Kp=p(1); %proportional gain
Ki=p(2); %integral gain
tau_max=p(3);

err = ucmd-u;
tau = Kp*err + Ki*interr;
%enforce limit
tau(tau > tau_max) = tau_max;
tau(tau < -tau_max) = -tau_max;

if nargout > 2
    dtau_du = -Kp*ones(length(u),1);
    dtau_du(abs(tau) == tau_max) = 0;
    dtau_du = diag(dtau_du);
end
