function [h_im,h_line] = plotPathOnAerial(x,y,lat,lon,zoom,border,h_axis)

    if nargin < 7
        h_axis = [];
        if nargin < 6
            border = [];
        end
    end
    
    if isempty(h_axis)
        set(figure,'name','path on aerial image')
        hold on
        xlabel('x (m)')
        ylabel('y (m)')
    else
        axes(h_axis)
    end
        
    clon = center(lon);
    clat = center(lat);
    cx = center(x);
    cy = center(y);
    
    url_fh = @googleStaticMapUrl;
%     url_fh = @bingStaticMapUrl;
    
    n = length(lat);
    I = floor(linspace(1,n,2)); %2 for just start & end
    
    url = feval(url_fh,lat(I),lon(I),zoom,clat,clon);
%     url = feval(url_fh,[],[],zoom,clat,clon); %no markers

    
    fn = '_autosave/staticmap.png';
    urlwrite(url,fn);
    [im,map]=imread(fn,'png');
    
    circ = 2*pi*6378137; %Earth circumference in meters
    metersPerPixel = (cos(clat * pi/180) * circ) / (256 * 2^zoom);
    
    %pixel coordinates of x=0,y=0
    x0 = size(im,2)/2 - cx/metersPerPixel;
    y0 = size(im,1)/2 + cy/metersPerPixel;
    
    xim = ((1:size(im,2)) - x0)*metersPerPixel;
    yim = ((size(im,1):-1:1) - (size(im,1)-y0))*metersPerPixel; %positive y is up
    
    set(gcf,'Colormap',map)
    h_im = image(xim, yim, im, 'CDataMapping', 'direct');
    set(gca,'YDir','normal')
    
    if nargout > 1
        h_line = plot(x,y);
    end
    
    if isempty(h_axis)
        axis equal
        makeLegible(14)
    end
    
    if ~isempty(border)
        %tight square axes
        d = max(abs([x-cx; y-cy]));
        d = d + border; %border in meters
        axis([[-d d]+cx [-d d]+cy]);
    end
    
    
end

function c = center(x)
    c = (min(x)+max(x))/2;
end


