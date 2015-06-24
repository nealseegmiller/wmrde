function [lat,lon] = UTMToLatLon(x,y,hemisphere,zone)

mstruct = getUTMZoneInfo(hemisphere,zone);
[lat,lon] = convertTransverseMercator(x,y,mstruct,true);