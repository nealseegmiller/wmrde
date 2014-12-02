function makeBoxHgt(hgt,lx,ly,lz,offset,C)
%create a box-shaped patch
%INPUTS
%hgt:       hgtransform object to parent to
%lx,ly,lz:  box dimensions
%offset:    3x1, translation of box center, can be empty [] for center at origin
%C:         color

if nargin < 6
    C = [1 1 0];
    if nargin < 5
        offset = [];
    end
end

vertices = [0  0  0;
            lx 0  0;
            lx ly 0;
            0  ly 0;
            0  0  lz;
            lx 0  lz;
            lx ly lz;
            0  ly lz];

vertices = vertices - repmat([lx ly lz]/2,8,1);

if ~isempty(offset)
    vertices = vertices + repmat(offset',size(vertices,1),1);
end

faces = [1 2 6 5;
         2 3 7 6;
         3 4 8 7;
         4 1 5 8;
         1 2 3 4;
         5 6 7 8];

h = patch('Vertices',vertices,'Faces',faces,'FaceColor',C);
set(h,'Parent',hgt)
% set(h,'FaceAlpha',.25)





            