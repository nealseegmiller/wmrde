function [h,surfidx] = surfacesHeight(surfaces,pts)
%INPUTS
%surfaces:      1 x ns cell array of surface objects
%pts:           2 x n, x,y point locations in world coords
%OUTPUTS
%h:             1 x n, height the highest surface at each point
%surfidx:       1 x n, index of surface used to compute h
%               h & surfidx are NaN if point is not within the boundary of any surface

n = size(pts,2);
h = -Inf*ones(1,n);
surfidx = NaN(1,n);

for i = 1:length(surfaces)
    h_ = surfaceHeight(surfaces{i},pts);
    I = find(h_ > h);
    h(I) = h_(I);
    surfidx(I) = i;
end
h(isnan(surfidx)) = NaN;