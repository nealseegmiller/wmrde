function im = loadTextureImage(filename,ppm)
%im, scaled image
%assumes square images (in pixels & meters)

%get image size L in meters
switch filename
    case 'cracked_dirt.jpg'
        L = 2.0;
    case 'dead_grass.jpg'
        L = 2.0;
    case 'desert_sand.jpg'
        L = 1.5;
    case 'wood_grain.jpg'
        L = 1.0;
    case 'green_grass.jpg'
        L = 1.0;
    otherwise
        
end

im = imread(filename);
scale=L/size(im,1)*ppm;
if scale > 1
    disp('WARNING, texture image resize scale > 1, recommended to decrease resolution')
end
im = imresize(im,scale);