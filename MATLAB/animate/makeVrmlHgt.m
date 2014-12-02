function makeVrmlHgt(hgt, filename, HT, flip_dim, draw_edges, fix_lighting, alpha)
%create a patch from a VRML 2.0/97 file
%INPUTS
%hgt:           hgtransform object to parent to
%filename:      VRML file name (.wrl)
%HT:            4x4 homogeneous transform to apply to vertices, can be empty []
%flip_dim:      0 for no flip, else {1,2,3} to flip vertices along x,y,z axis
%draw_edges:    if true, draw lines on edges where normals of two adjoining faces differ sufficiently
%fix_lighting:  if true, duplicate vertices for multiple normal directions
%alpha:         transparency, 0 to 1 (clear to opaque). if 0 no triangles are drawn


if nargin < 7
    alpha=1;
    if nargin < 6
        fix_lighting = true;
        if nargin < 5
            draw_edges = false;
            if nargin < 4
                flip_dim = 0;
                if nargin < 3
                    HT = [];
                end
            end
        end
    end
end



[vrml.coord, vrml.coord_idx, vrml.color, vrml.color_idx, vrml.normal, vrml.normal_idx] = parseVrml(filename);

np = length(vrml.coord); %number of patches

for pno = 1:np %part number

    nt = size(vrml.coord_idx{pno},1); %number of triangles
    
    vertices = vrml.coord{pno};
    faces = vrml.coord_idx{pno};
    
    if ~isempty(HT)
        vertices = applyHT(HT,vertices')';
    end

    if flip_dim ~= 0
        vertices(:,flip_dim) = -vertices(:,flip_dim);

        %no longer right-handed, flip cw/ccw
        %is this necessary? affects lighting?
%         faces(:,[2 3]) = faces(:,[3 2]);
    end
    
    if alpha > 0

        colors = vrml.color{pno}(vrml.color_idx{pno},:); %per face colors

        if fix_lighting
            tmp = vrml.coord_idx{pno}';
            vertices_ = vertices(tmp(:),:);
            faces_ = reshape(1:(nt*3),3,nt)';

            %not necessary to set VertexNormals, MATLAB computes them automatically
%             tmp = ones(3,1)*(1:nt);
%             normals = vrml.normal{pno}(vrml.normal_idx{pno}(tmp(:)),:);      
        else        
            vertices_ = vertices;
            faces_ = faces;
        end
    
        h = patch('Parent',hgt,'Vertices',vertices_,'Faces',faces_,'FaceVertexCData',colors);
        set(h,'EdgeColor','none','FaceColor','flat','FaceLighting','flat','FaceAlpha',alpha);
        
%         if fixlighting
%             set(h,'VertexNormals',-normals);
%         end

        if 0
            %DEBUGGING, display vertex normals
            N = -get(h,'VertexNormals')'; %3 x n
            N = N./( ones(3,1)*sqrt(sum(N.^2,1)) ); %normalize

            hq = drawQuivers(vertices_',N/20);
            set(hq,'Color','red')
        end
    end
    

    

    if draw_edges
        tol = 15*pi/180; %deg
        %DRAW LINES
        %every edge is shared by two triangle faces
        %only draw edges for which the normals of those two triangles differ sufficiently
        
        nv = size(vrml.coord{pno},1);
        
        normal_idx_1 = zeros(nv,nv); %normal index for 1st triangle that contains edge
        normal_idx_2 = zeros(nv,nv); %2nd

        efaces = []; %nx2, edge faces

        for tno = 1:nt %loop over triangles
            %loop over triangle edges
            for vno_1 = 1:3 %1st vertex no (in the triangle)
                vno_2 = vno_1 + 1; %2nd vertex no
                if (vno_2 > 3)
                    vno_2 = rem(vno_2,3); 
                end
                
                %vertex indices in the entire list
                v1 = vrml.coord_idx{pno}(tno,vno_1);
                v2 = vrml.coord_idx{pno}(tno,vno_2);
                
                %sort
                tmp = sort([v1 v2]);
                v1 = tmp(1); 
                v2 = tmp(2);
                
                if normal_idx_1(v1,v2) == 0
                    normal_idx_1(v1,v2) = vrml.normal_idx{pno}(tno);
                else
                    normal_idx_2(v1,v2) = vrml.normal_idx{pno}(tno);

                    %check angle between normal vectors, draw line if > tol
                    N1 = vrml.normal{pno}(normal_idx_1(v1,v2),:)';
                    N2 = vrml.normal{pno}(normal_idx_2(v1,v2),:)';
                    angle = acos(N1'*N2);
                    if abs(angle) > tol
                        efaces = [efaces; v1 v2]; %#ok<AGROW>
                    end 
                end
            end
        end

        patch('Parent',hgt,'Vertices',vertices,'Faces',efaces,'FaceColor','none');

    end
    


end




