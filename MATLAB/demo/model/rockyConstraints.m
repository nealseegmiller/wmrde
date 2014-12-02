function [c, Jc, f, df_djd, df_djr] = rockyConstraints(mdl, jd, jr)
%Rocky differential averaging constraint

nj = length(jd); %num. joints
% diff_fi = namesToInds(mdl,{'D1','D2'}); %differential frame indices
diff_fi = [2 8];
diff_ji = diff_fi - 1; %joint indices

c = sum(jd(diff_ji)); %angles should be equal magnitude, opposite sign

Jc = zeros(1,nj);
Jc(diff_ji) = 1;

if nargin > 2
    Kp = 2250;
    Kd = Kp/20;
    cd = sum(jr(diff_ji)); %d/dt c
    f = -(Kp*c + Kd*cd);

    df_djd = -Kp*Jc;
    df_djr = -Kd*Jc; %dcd/djr = dc/djd

end








