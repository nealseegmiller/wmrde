function [northing,easting,lonzone] = LatLonToUTM(lat,lon)


%datum = WGS 84
a = 6378137;
b = 6356752.3142;

k0 = .9996;

e = sqrt(1-(b/a)^2);
eisq = e*e/(1-e*e);

n=(a-b)/(a+b);
A0=a*(1-n+(5*n*n/4)*(1-n) +(81*n^4/64)*(1-n));
B0=(3*a*n/2)*(1 - n - (7*n*n/8)*(1-n) + 55*n^4/64);
C0=(15*a*n*n/16)*(1 - n +(3*n*n/4)*(1-n));
D0=(35*a*n^3/48)*(1 - n + 11*n*n/16);
E0=(315*a*n^4/51)*(1-n);


lonzone = 31+floor(lon/6);
lonzonecm = 6*lonzone-183;
deltalon = (lon-lonzonecm)*pi/180; %rad
lat = lat*pi/180; %convert to radians

% rho = a*(1-e*e)/((1-(e*sin(lat))^2)^(3/2)); %r curv 1
nu = a./((1-(e*sin(lat)).^2).^(1/2)); %r curv 2
MeridionalArc = A0*lat - B0*sin(2*lat) + C0*sin(4*lat) - D0*sin(6*lat) + E0*sin(8*lat);
ki = MeridionalArc*k0;
kii = nu.*sin(lat).*cos(lat)/2;
kiii = ((nu.*sin(lat).*cos(lat).^3)/24).*(5-tan(lat).^2+9*eisq*cos(lat).^2+4*eisq^2*cos(lat).^4)*k0;
kiv = nu.*cos(lat)*k0;
kv = (cos(lat)).^3.*(nu/6).*(1-tan(lat).^2+eisq*cos(lat).^2)*k0;
% A6 = ((deltalon)^6*nu*sin(lat)*cos(lat)^5/720)*(61-58*tan(lat)^2+tan(lat)^4+270*eisq*cos(lat)^2-330*eisq*sin(lat)^2)*k0;
northing = (ki+kii.*deltalon.*deltalon+kiii.*deltalon.^4);

northing(northing < 0) = 1e7+northing(northing < 0);

easting = 500000+(kiv.*deltalon+kv.*deltalon.^3);

