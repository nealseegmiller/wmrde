function p = wgcParams(fh,Kp)
%set parameters for the wheel-ground contact model
%INPUTS
%fh, wheel-ground contact model function handle
%Kp, tire stiffness (total mass)*gravity/((number of wheels)*compression)
%OUTPUTS
%p, parameter vector

Kd = Kp/20;

%rolling resistance
% Crr=-.05;
Crr=0; %DEBUGGING



switch func2str(fh)
    case 'odeWgc'
        s = 1e-3; %N to kN
       
        mu = 1.0;
%         mu = Inf;
        
        C = -.1*Kp;
        
        p = [Kp*s Kd*s mu mu C*s C*s Crr*s]';
        
    case 'pacejkaWgc'
        mu = 1.0;
%         mu = 0.6;
        
        %from Brach & Brach, {B,C,D,E,K}
        clon = [1/15 1.5 1 .3 100];
        clat = [8/75 1.5 1 .6 100];
        
%         clon(4)=1; clat(4)=1; %DEBUGGING, E=1 no instability
        
        p = [Kp Kd mu, clon, clat Crr]';
    case {'ishigamiWgc','ishigamiLUTWgc'}
        
        p = [Kd, 0.18/2, 0.11, 0.8, 37.2*pi/180, 26.4*pi/180, 1.37e3, 8.14e5, 1.0, 0.4, 0.15, 1.6e3*9.81, 1];
%         p = [p, .043, .036, .020, .013]'; %Ishigami
        p = [p, .025]'; %Jia, isotropic
        
    otherwise
        disp([mfilename ': wheel-ground contact model function handle not recognized'])
end


