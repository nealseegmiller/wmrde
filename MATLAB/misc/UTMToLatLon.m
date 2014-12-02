function [lat,lon] = UTMToLatLon(northing,easting,hemisphere,lonzone)
%northing, nx1 vector
%easting, nx1 vector
%hemisphere, 1xn string 'N','S'
%lonzone, nx1 vector (integers)

% http://www.uwgb.edu/dutchs/usefuldata/utmformulas.htm#Spreadsheet
% http://www.uwgb.edu/dutchs/usefuldata/UTMConversions1.xls
  
%UTM longitudinal zone map:
% http://www.dmap.co.uk/utmworld.htm

%datum = WGS 84
a = 6378137;
b = 6356752.3142;

k0 = .9996;

e = sqrt(1-(b/a)^2);

ei = (1-(1-e*e)^(1/2))/(1+(1-e*e)^(1/2));
c1 = 3*ei/2-27*ei^3/32;
c2 = 21*ei^2/16-55*ei^4/32;
c3 = 151*ei^3/96;
c4 = 1097*ei^4/512;

eisq = e*e/(1-e*e);


east_prime = 5e5-easting;
arc_length = northing/k0;
mu = arc_length/(a*(1-e^2/4-3*e^4/64-5*e^6/256));
phi=mu+c1*sin(2*mu)+c2*sin(4*mu)+c3*sin(6*mu)+c4*sin(8*mu); %footprint latitude

C1 = eisq*cos(phi).^2;
T1 = tan(phi).^2;
N1 = a./(1-(e*sin(phi)).^2).^(1/2);
R1 = a*(1-e*e)./(1-(e*sin(phi)).^2).^(3/2);

D = east_prime./(N1*k0);
Fact1 = N1.*tan(phi)./R1;
Fact2 = D.*D/2;
Fact3 = (5+3*T1+10*C1-4*C1.*C1-9*eisq).*D.^4/24;
Fact4 = (61+90*T1+298*C1+45*T1.*T1-252*eisq-3*C1.*C1).*D.^6/720;
LoFact1 = D;
LoFact2 = (1+2*T1+C1).*D.^3/6;
LoFact3 = (5-2*C1+28*T1-3*C1.^2+8*eisq+24*T1.^2).*D.^5/120;
DeltaLon = (LoFact1-LoFact2+LoFact3)./cos(phi);
ZoneCM = 6*lonzone-183;
lat = 180*(phi-Fact1.*(Fact2+Fact3+Fact4))./pi;

issouth= hemisphere=='S';
lat(issouth)=-lat(issouth);

lon = ZoneCM-DeltaLon*180/pi;







