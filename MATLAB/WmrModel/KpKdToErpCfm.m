function [erp,cfm] = KpKdToErpCfm(kp,kd,dt)

%conversion per ODE documentation
%http://ode-wiki.org/wiki/index.php?title=Manual:_Concepts

erp = dt*kp/(dt*kp + kd);
cfm = 1/(dt*kp + kd);