function url = googleStaticMapUrl(lat,lon,zoom,clat,clon)

% https://developers.google.com/maps/documentation/staticmaps/
% alternative:
% http://earthexplorer.usgs.gov/


%get url for google static maps api
fmt='%.6f'; %format for lat,lon coords

url = ['http://maps.google.com/maps/api/staticmap?center=' num2str(clat,fmt) ',' num2str(clon,fmt) ...
       '&zoom=',num2str(zoom),'&size=640x640&maptype=satellite'];

%sizes: {tiny, mid, small}
%colors: {black, brown, green, purple, yellow, blue, gray, orange, red, white}
n = length(lat);
for i=1:n
    if i==1, color = 'green'; 
    elseif i==n, color = 'red';
    else color = 'black';
    end
    url = [url, '&markers=size:tiny%7Ccolor:' color '%7C' num2str(lat(i),fmt) ',' num2str(lon(i),fmt)]; %#ok<AGROW>
end

url = [url, '&sensor=false'];

