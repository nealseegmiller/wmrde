function mstruct = getUTMZoneInfo(in1,in2)
% INPUTS
% in1,in2 = hemisphere,zone OR lat,lon
% 
% source:
% http://en.wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system#Simplified_formulas
% map unit = meters
% datum = WGS 84

if ischar(in1)
    hemisphere = in1; 
    zone = in2;
else
    lat = in1;
    lon = in2;
    
    if lat > 0
        hemisphere = 'N';
    else
        hemisphere = 'S';
    end
    zone = 31+floor(lon/6);
end
%copy to struct
mstruct.hemisphere = hemisphere;
mstruct.zone = zone;

a = 6378137; %equatorial radius
b = 6356752.3142; %polar radius
ec = sqrt(1-(b/a)^2); %eccentricity

mstruct.geoid = [a ec];

if hemisphere == 'N'
    mstruct.falsenorthing = 0;
elseif hemisphere == 'S'
    mstruct.falsenorthing = 1e7;
end
mstruct.falseeasting = 5e5;
mstruct.scalefactor = .9996;

ZoneCM = 6*zone-183;
mstruct.origin = [0 ZoneCM 0];


