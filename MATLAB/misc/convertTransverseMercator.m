function [out1,out2] = convertTransverseMercator(in1,in2,mstruct,inverse)
% Convert between Transverse Mercator coordinates (i.e. UTM or State Plane) and latitude/longitude
% INPUT
% in1:      n x 1
% in2:      n x 1
% mstruct:  struct with zone information (see below)
% inverse:  flag
%   if false do forward conversion from lat,lon to x,y (easting,northing)
%   if true do inverse conversion from x,y to lat,lon
% OUTPUT
% out1:     n x 1
% out2:     n x 1
% 
% source of formulas:
% http://www.uwgb.edu/dutchs/usefuldata/utmformulas.htm
% alternative source:
% http://pubs.er.usgs.gov/publication/pp1395, p 60-63
% to validate output:
% http://www.rcn.montana.edu/resources/converter.aspx
% http://www.earthpoint.us/StatePlane.aspx

angRad = pi/180; %multiply by this to convert deg to rad

%parse mstruct
a = mstruct.geoid(1);                   %equatorial radius, map units
e = mstruct.geoid(2);                   %eccentricity of ellipsoid
East0 = mstruct.falseeasting;           %false easting of central meridian, map units
North0 = mstruct.falsenorthing;         %false northing, map units
k0 = mstruct.scalefactor;
lat0 = mstruct.origin(1) * angRad;      %DEBUGGING
ZoneCM = mstruct.origin(2) * angRad;    %central meridian


ei2 = e*e/(1-e*e); %e prime squared
    
if inverse
    %from x,y to lat,lon

    x = in1;
    y = in2;
    
    %compute the footprint latitude
    M = (y-North0)/k0; %Meridional Arc
    
    M = M + calcMeridionalArc(lat0,a,e);
    
    mu = M/(a*(1 - e^2/4 - 3*e^4/64 - 5*e^6/256));
    e1 = (1 - sqrt(1-e^2))/(1 + sqrt(1-e^2));
    
    J1 = 3*e1/2 - 27*e1^3/32;
    J2 = 21*e1^2/16 - 55*e1^4/32;
    J3 = 151*e1^3/96;
    J4 = 1097*e1^4/512;

    fp = mu + J1*sin(2*mu) + J2*sin(4*mu) + J3*sin(6*mu) + J4*sin(8*mu);

    %compute lat,lon
    C1 = ei2*cos(fp).^2;
    T1 = tan(fp).^2;
    R1 = a*(1-e^2) ./ (1-(e*sin(fp)).^2).^(3/2);
    N1 = a./sqrt(1-(e*sin(fp)).^2); %c.f. nu in forward conversion
    D = (x - East0)./(N1*k0);
    
    Q1 = N1.*tan(fp)./R1;
    Q2 = D.*D/2;
    Q3 = (5 + 3*T1 + 10*C1 - 4*C1.*C1 - 9*ei2).*(D.^4)/24;
    Q4 = (61 + 90*T1 + 298*C1 + 45*T1.*T1 - 3*C1.*C1 - 252*ei2 ).*(D.^6)/720;
    
    lat = fp - Q1.*(Q2-Q3+Q4);
    
    Q5 = D;
    Q6 = (1 + 2*T1 + C1).*(D.^3)/6;
    Q7 = (5 - 2*C1 + 28*T1 - 3*C1.^2 + 8*ei2 + 24*T1.^2).*(D.^5)/120;
    
    lon = ZoneCM + (Q5 - Q6 + Q7)./cos(fp);
    
    out1 = lat / angRad; %convert to deg;
    out2 = lon / angRad;
    
else %forward

    %from lat,lon to x,y
   
    lat = in1 * angRad; %convert to rad
    lon = in2 * angRad;

    M = calcMeridionalArc(lat,a,e) - calcMeridionalArc(lat0,a,e);

    %radius of curvature of the earth perpendicular to the meridian plane
    nu = a./sqrt((1-(e*sin(lat)).^2)); 
    
    p = lon - ZoneCM;
    p2 = p.*p;
    p3 = p2.*p;
    p4 = p3.*p;
    
    K1 = M*k0;
%     K2 = k0*nu.*sin(lat).*cos(lat)/2;
    K2 = k0*nu.*sin(2*lat)/4; %equivalent
    K3 = (k0*nu.*sin(lat).*cos(lat).^3/24) .* (5 - tan(lat).^2 + 9*ei2*cos(lat).^2 + 4*ei2^2*cos(lat).^4);
    
    y = K1 + K2.*p2 + K3.*p4 + North0;
    
    K4 = k0*nu.*cos(lat);
    K5 = (k0*nu.*cos(lat).^3/6) .* (1 - tan(lat).^2 + ei2*cos(lat).^2);

    x = K4.*p + K5.*p3 + East0;
    
    out1 = x;
    out2 = y;
    
end

end


function M = calcMeridionalArc(phi,a,e)
    %compute the Meridional Arc
    %http://en.wikipedia.org/wiki/Meridian_arc
    
    persistent A0 B0 C0 D0 E0
    if isempty(A0)
        if 1
            
            b = a*sqrt(1-e^2); %polar radius
            n = (a-b)/(a+b);
            n2 = n*n;
            n3 = n2*n;
            n4 = n3*n;
            n5 = n4*n;

            A0 = a*(1 - n + (5/4)*(n2-n3) + (81/64)*(n4-n5));
            B0 = (3/2*a*n)*(1 - n + (7/8)*(n2-n3) + (55/54)*(n4-n5));
            C0 = (15/16*a*n2)*(1 - n + (3/4)*(n2-n3));
            D0 = (35/48*a*n3)*(1 - n + (11/16)*(n2-n3));
            E0 = (315/512*a*n4)*(1-n);

        else

            %USGS
            e2 = e*e;
            e4 = e2*e2;
            e6 = e4*e2;

            A0 = a*(1 - e2/4 - 3*e4/64 - 5*e6/256);
            B0 = a*(3*e2/8 + 3*e4/32 + 45*e6/1024);
            C0 = a*(15*e4/256 + 45*e6/1024);
            D0 = a*35*e6/3072;
            E0 = 0;

        end
    end
    
    M = A0*phi - B0*sin(2*phi) + C0*sin(4*phi) - D0*sin(6*phi) + E0*sin(8*phi);
end
