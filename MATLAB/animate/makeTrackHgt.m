function makeTrackHgt(hgt,rad1,rad2,L,width,offset,enc,C)

if nargin < 8
    C = [1 1 0];
    if nargin < 7
        enc = 1;
        if nargin < 6
            offset = [];
        end
    end
end


%make a circle for the first sprocket
circle_nv = 18;

circ = pointsOnCircle(rad1,circle_nv,1);
poly_pts = circ;

%make a trapezoid
trap = [0 rad1; L rad2; L -rad2; 0 -rad1; 0 rad1]';

[xu,yu] = polybool('union', poly_pts(1,:), poly_pts(2,:), trap(1,:), trap(2,:));
poly_pts = [xu; yu];

%make a circle for the 2nd sprocket
circ = pointsOnCircle(rad2,circle_nv,1);
circ(1,:) = circ(1,:) + L;

[xu,yu] = polybool('union', poly_pts(1,:), poly_pts(2,:), circ(1,:), circ(2,:));
poly_pts = [xu; yu];


%cull vertices that are too close to neighboring vertices
arc_length = 2*pi*rad1/circle_nv;
tol = arc_length/4;

j = 2;
while j < size(poly_pts,2)
    if norm(poly_pts(:,j) - poly_pts(:,j-1)) < tol
        poly_pts(:,j) = [];
    else
        j=j+1;
    end
end


%extrude the polygon
makePolyExtrudeHgt(hgt,poly_pts,width,2,offset,enc,false,C);




