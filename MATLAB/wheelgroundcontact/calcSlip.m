function [s,alpha,ds_dvx,ds_dRw,dalpha_dvx,dalpha_dvy,dalpha_dRw] = calcSlip(vx,vy,Rw,method)
%INPUTS
%vx lon velocity of contact point
%vy lat velocity of contact point
%Rw = R*omega, wheel radius * wheel angular velocity
%method (see below)
%OUTPUTS
%s, slip ratio
%alpha, slip angle
%computed according to method:
%method=1; (Wong eq 1.6, SAE, Brach & Brach) s = (Rw-Vx)/Vx = -vx/(vx + Rw)
%method=2; (Wong eq 1.5, Iagnemma) s = (Rw-Vx)/Rw = -vx/Rw
%(all methods) alpha = atan(vy/Vx)

%Ishigami, Jia propose switching between methods 1 & 2 based on braking/driving
%but the shear displacement equations they use for jx,jy assume method 2

Vx = vx + Rw; %lon velocity of wheel center

%slip angle not well defined when Vx near zero
%TODO, a better way? derivatives okay?
tol=1e-3; %(m/s)
Vx(Vx<tol & Vx>=0) = tol;
Vx(Vx>-tol & Vx<0) = -tol;

if method==1
    s = -vx./Vx;
elseif method==2
    s = -vx./Rw;
end

%slip angle
tan_alpha = vy./Vx;
alpha = atan(tan_alpha);

if nargout > 2
    %compute derivatives
    
    Vx2 = Vx.^2; %precompute
    
    if method==1
        ds_dvx = -1./Vx + vx./Vx2;
        ds_dRw = vx./Vx2;
        
    elseif method==2
        ds_dvx = -1./Rw;
        ds_dRw = vx./Rw.^2;
        
    end

    temp1 = 1./(1+(tan_alpha).^2); %d/dx atan(x) = 1/(1+x^2)
    
    dalpha_dvx = temp1.*(-vy./Vx2);
    dalpha_dvy = temp1./Vx;
    dalpha_dRw = dalpha_dvx;

end





