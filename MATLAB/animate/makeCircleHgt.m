function makeCircleHgt(hgt,rad,ax,offset,nv,C)

if nargin < 6
    C = [0 1 1]; %default color
    if nargin < 5
        nv = 18;
        if nargin < 4
            offset = [];
        end
    end
end

vertices = pointsOnCircle(rad,nv,1,ax)'; 
np = nv+1;

if ~isempty(offset)
    vertices = vertices + repmat(offset',size(vertices,1),1);
end

h = patch('Vertices',vertices,'Faces',1:np,'FaceColor',C);
set(h,'Parent',hgt)

%make an line so rotation is visible
h = line([offset(1) vertices(1,1)],[offset(2) vertices(1,2)],[offset(3) vertices(1,3)],'Color',[0 0 0]);
set(h,'Parent',hgt)