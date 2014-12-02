function N = surfacesNormal(surfaces,pts,surfidx)
%INPUTS
%pts:       2 x n, x,y point locations in world coords
%surfidx:   1 x n, index of surface in surfaces to use
%OUTPUTS
%N:         3 x n, surface normals at each point

n = size(pts,2);
N = NaN*ones(3,n);

for i = 1:length(surfaces)
    N(:,surfidx==i) = surfaceNormal(surfaces{i},pts(:,surfidx==i));
end