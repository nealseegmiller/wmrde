function [c, Jc, f, df_djd, df_djr] = recbotConstraints(mdl, jd, jr)
%Gator suspension constraints

nj = length(jd); %num. joints
% fi = names_to_inds(obj,{'suspL','suspR'}); %suspension frame indices
fi = 2:3;
I = fi - 1; %suspension joint indices

%TODO, kinematic sim

c = jd(I);

Jc = eye(nj);
Jc = Jc(I,:); %dc/djd

if nargin > 2
    Kp = 3e4;
    Kd = Kp/20;
    cd = jr(I); %d/dt c
    f = -(Kp*c + Kd*cd);

    df_djd = -Kp*Jc;
    df_djr = -Kd*Jc; %dcd/djr = dc/djd

end

