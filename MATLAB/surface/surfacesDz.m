function [dz,surfidx] = surfacesDz(surfaces,pts)
%INPUTS
%surfaces:  1 x ns cell array of surface objects
%pts:       3 x n, x,y,z point locations in world coords
%OUTPUTS
%dz:        1 x n, contact height errors for each point. Distance from point to surface along surface normal. Use surface for which dz is smallest.
%surfidx:   1 x n, index of surface used to compute dz
%           dz & surfidx are NaN if point is not within the boundary of any surface

n = size(pts,2);
dz = Inf*ones(1,n);
surfidx = NaN(1,n);

for i = 1:length(surfaces)
    dz_ = surfaceDz(surfaces{i},pts);
    I = find(dz_ < dz);
    dz(I) = dz_(I);
    surfidx(I) = i;
end
dz(isnan(surfidx)) = NaN;
