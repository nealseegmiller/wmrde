function makeLineExtrudeHgt(hgt,line_pts,thk,height,ax,offset,enc,C)
%create an extruded polygon shaped patch derived from a 2D line
%fails if self intersecting?
%INPUTS
%hgt:       hgtransform object to parent to
%line_pts:  2xn matrix of points on line
%thk:       line thickness
%height:    extrusion height
%ax:        extrusion axis {1,2,3}
%offset:    3x1, translation of center, can be empty [] for center at origin
%enc:       if true, enclose with faces
%C:         color

if nargin < 8
    C = [1 1 0];
    if nargin < 7
        enc = 1;
        if nargin < 6
            offset = [];
        end
    end
end

opts.debug_plot = 0;

%make a circle polygon around the first point
i = 1;
rad = thk/2;
circle_nv = 18;

circ_ = pointsOnCircle(rad,circle_nv,1);
poly_pts(1,:) = circ_(1,:) + line_pts(1,i); %move circle
poly_pts(2,:) = circ_(2,:) + line_pts(2,i);

if opts.debug_plot
    %DEBUGGING
    figure;
    axis equal; hold on
    colors = 'rgbkmc';
end

for i = 2:size(line_pts,2)
    pt1 = line_pts(:,i-1);
    pt2 = line_pts(:,i);
    
    %make a rectangle polygon from the edge between pts i-1 and i
    %edge normal unit vector
    normal = [0 -1; 1 0]*(pt2-pt1)/norm(pt2-pt1);
    nr = normal*rad;
    rect = [pt1 + nr, pt2 + nr, pt2 - nr, pt1 - nr, pt1 + nr];
    
    %union with vertices
    [xu,yu] = polybool('union', poly_pts(1,:), poly_pts(2,:), rect(1,:), rect(2,:));
    poly_pts = [xu; yu];
    
    %make a circle polygon around point i
    circ = circ_;
    circ(1,:) = circ(1,:) + pt2(1);
    circ(2,:) = circ(2,:) + pt2(2);
    
    %union with vertices
    [xu,yu] = polybool('union', poly_pts(1,:), poly_pts(2,:), circ(1,:), circ(2,:));
    poly_pts = [xu; yu];
    
    if opts.debug_plot
        plot(poly_pts(1,:),poly_pts(2,:),[colors(i) '.-']); %DEBUGGING
    end
    
end

%cull vertices that are too close to neighboring vertices
arc_length = 2*pi*rad/circle_nv;
tol = arc_length/4;
j = 2;
while j < size(poly_pts,2)
    if norm(poly_pts(:,j) - poly_pts(:,j-1)) < tol
        poly_pts(:,j) = [];
    else
        j=j+1;
    end
end

if opts.debug_plot
    plot(poly_pts(1,:),poly_pts(2,:), [colors(i+1) ':']) %DEBUGGING
end

%extrude the polygon
makePolyExtrudeHgt(hgt,poly_pts,height,ax,offset,enc,false,C);




