function url = bingStaticMapUrl(lat,lon,zoom,clat,clon)

% http://msdn.microsoft.com/en-us/library/ff701724.aspx
% http://msdn.microsoft.com/en-us/library/bb259689.aspx 
% ground resolution (meters/pixel) = (cos(latitude * pi/180) * 2 * pi * 6378137 meters) / (256 * 2^level pixels)


%get url for bing static maps api
fmt='%.6f'; %format for lat,lon coords

url = ['http://dev.virtualearth.net/REST/v1/Imagery/Map/Aerial/' num2str(clat,fmt) ',' num2str(clon,fmt) '/' num2str(zoom) '?mapSize=800,800' ...
       '&format=png' ...
       '&key=' bingmapskey()];

for i=1:length(lat)
    url = [url, '&pushpin=' num2str(lat(i),fmt) ',' num2str(lon(i),fmt) ';9']; %#ok<AGROW>
end
