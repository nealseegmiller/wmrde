%symbolic calculation of derivatives for 
% Pacejka Nicolas-Comstock wheel-ground contact model

clear_all
close all

syms('pi')
syms('mu','Bs','Cs','Ds','Es','Ks','Ba','Ca','Da','Ea','Ka') %parameters
syms('s','alpha','fz') %inputs
syms('sgn')

%Pacejka, forces for independent lon/lat slip
Fxs = Ds*sin(Cs*atan(Bs*(1-Es)*Ks*s + Es*atan(Bs*Ks*s)));
temp1 = (2/pi)*alpha;
Fya = Da*sin(Ca*atan(Ba*(1-Ea)*Ka*temp1 + Ea*atan(Ba*Ka*temp1)));

%Nicolas Comstock, forces combined lon/lat slip
s2 = s.*s;
Fxs2 = Fxs.*Fxs;
Fya2 = Fya.*Fya;
tan_alpha = tan(alpha);
tan_alpha2 = tan_alpha.*tan_alpha;
temp2 = abs(Fxs.*Fya) ./ sqrt( s2.*Fya2 + Fxs2.*tan_alpha2 );

fx = sgn .* temp2 .* s .* (mu*fz);
fy = -sgn .* temp2 .* tan_alpha .* (mu*fz);

dfx_ds = diff(fx,s);
dfy_ds = diff(fy,s);
dfx_dalpha = diff(fx,alpha);
dfy_dalpha = diff(fy,alpha);

matlabFunction(dfx_ds, dfy_ds, dfx_dalpha, dfy_dalpha, 'file', 'pacejkaDerivatives.m')
ccode([dfx_ds dfy_ds, dfx_dalpha, dfy_dalpha], 'file', 'pacejkaDerivatives.cpp')

% matlabFunction(limit(temp2,s,0,'right'))
% matlabFunction(limit(temp2,alpha,0,'right'))

