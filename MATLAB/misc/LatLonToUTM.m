function [x,y] = LatLonToUTM(lat,lon)
%assumes all points are in the same zone!

mstruct = getUTMZoneInfo(lat(1),lon(1));
[x,y] = convertTransverseMercator(lat,lon,mstruct,false);
