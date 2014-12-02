function makeCylinderHgt(hgt,rad,height,ax,offset,nv,enc,C)
%create a cyclinder shaped patch
%INPUTS
%hgt:       hgtransform object to parent to
%rad:       cylinder radius
%height:    cylinder height
%ax:        {1,2,3} axis of symmetry
%offset:    3x1, translation of cylinder center, can be empty [] for center at origin
%nv:        number of vertices in circle
%enc:       if true, enclose cylinder with circular faces
%C:         color

if nargin < 8
    C = [0 1 1]; %default color
    if nargin < 7
        enc = 1; %enclose
        if nargin < 6
            nv = 18; %number of vertices on circle
            if nargin < 5
                offset = [];
            end
        end
    end
end


pts = pointsOnCircle(rad,nv,1); %output is 2x(nv+1) because closed polygon
np = nv+1;

vertices = makePolyExtrudeHgt(hgt,pts,height,ax,offset,enc,0,C);

%make a line so rotation is visible
p = vertices([1, np+1],:);
if enc
    p = [offset'; p; offset'];
    p(1,ax) = p(1,ax) - height/2;
    p(4,ax) = p(4,ax) + height/2;
end
h = line(p(:,1),p(:,2),p(:,3),'Color',[0 0 0]);
set(h,'Parent',hgt);

