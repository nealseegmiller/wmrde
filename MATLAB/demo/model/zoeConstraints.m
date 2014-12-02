function [c, Jc, f, df_djd, df_djr] = zoeConstraints(mdl, jd, jr)
%Zoe axle roll averaging constraint

nj = length(jd); %num. joints
% fi = names_to_inds(mdl,{'r1','r2'}); %axle roll frame indices
roll_fi = [2 6];
roll_ji = roll_fi - 1; %differential joint indices

c = sum(jd(roll_ji)); %angles should be equal magnitude, opposite sign

Jc = zeros(1,nj);
Jc(roll_ji) = 1;

if nargout > 2
    Kp = 22500;
    Kd = Kp/20;
    cd = sum(jr(roll_ji)); %d/dt c
    f = -(Kp*c + Kd*cd);

    df_djd = -Kp*Jc;
    df_djr = -Kd*Jc; %dcd/djr = dc/djd
end



