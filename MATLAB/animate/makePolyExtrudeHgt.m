function vertices = makePolyExtrudeHgt(hgt,pts,height,ax,offset,enc,show_edges,C)

%pts, 2xnp matrix specifying vertices of a closed polygon
%axis (1,2,3)
%   if axis == 1; pts = [Y; Z], extrude along x
%   if axis == 2; pts = [X; Z], extrude along y
%   if axis == 3; pts = [Y; Z], extrude along z

if nargin < 8
    C = [1 1 0]; %default color
    if nargin < 7
        show_edges = 1; %default, extruded edges shown by lines
        %top and bottom face edge always visible
        if nargin < 6
            enc = 1; %enclose the patch
            if nargin < 5
                offset = [];
            end
        end
    end
end
%default origin is midpoint in extrusion
offset(ax) = offset(ax) - height/2;
    
np = size(pts,2); %number of pts
I = 1:3;
I(ax) = [];

v1 = zeros(np,3); %vertices of bottom face
v1(:,I) = pts';
v2 = v1; %vertices of top face
v2(:,ax) = height;

%add the offset
if ~isempty(offset)
    v1 = v1 + repmat(offset',size(v1,1),1);
    v2 = v2 + repmat(offset',size(v2,1),1);
end

vertices = [v1; v2];

%make patch of extruded faces
faces = zeros(np-1,4);
for i = 1:(np-1)
    faces(i,:) = [i i+1 i+np+1 i+np];
end

if show_edges
    edgeC = [0 0 0];
else
    edgeC = 'none';
end

alpha=1; %TODO, make input

h = patch('Vertices',vertices,'Faces',faces,'FaceColor',C,'EdgeColor',edgeC,'FaceAlpha',alpha);
set(h,'Parent',hgt)


%make top and bottom faces
if enc
    C_ = C;
else
    C_ = 'none'; %still make patch just for edge
end

h = patch('Vertices',v1,'Faces',1:(np-1),'FaceColor',C_,'FaceAlpha',alpha);
set(h,'Parent',hgt)

h = patch('Vertices',v2,'Faces',1:(np-1),'FaceColor',C_,'FaceAlpha',alpha);
set(h,'Parent',hgt)
