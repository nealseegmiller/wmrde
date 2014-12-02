function [kp,kd] = ErpCfmToKpKd(erp,cfm,dt)

%DEBUGGING, check this!
kd = -(erp-1)/cfm;
kp = erp/(cfm*dt);

