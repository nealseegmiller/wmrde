classdef TriMeshSurf < Surface
    
    properties (GetAccess = 'public', SetAccess = 'private')
        
        vertices    % 3 x nv, in world coords
        indices     % nt x 3, vertex indices for triangles
        pec         % nt x 4, plane equation coefficients for each triangle, 
                    %pec * [x y z 1] = distance of point from surface, pec(1:3) is surface normal
        %inequality matrices, for determining if point is in bounds of triangle
        %one for each edge of triangle (vertex 1-2, 2-3, and 3-1)
        M12     % nt x 3, 
        M23     % nt x 3
        M31     % nt x 3

        %bounding box
        xmin
        xmax
        ymin
        ymax
        
    end
    
    properties (Dependent)
        nv          % number of vertices
        nt          % number of triangles
    end
    
    
    methods
        function obj = TriMeshSurf()
            obj.vertices = zeros(3,0);
            obj.indices = zeros(0,3);
            
            obj.M12 = zeros(0,3);
            obj.M23 = zeros(0,3);
            obj.M31 = zeros(0,3);
            
        end
        
        %Dependent access methods
        function out = get.nv(obj)
            out = size(obj.vertices,2);
        end
        function out = get.nt(obj)
            out = size(obj.indices,1);
        end
        
        function addVertices(obj, V)
            %V:     3 x n
            obj.vertices = [obj.vertices, V];
            
            %bounding box
            obj.xmin = min(obj.vertices(1,:));
            obj.xmax = max(obj.vertices(1,:));
            
            obj.ymin = min(obj.vertices(2,:));
            obj.ymax = max(obj.vertices(2,:));
        end
        
        %indices must be specified in clockwise order, viewed from above.
        %z component of normal must be > 0
        function addTriangle(obj, I)
            %I:     size 3, indices
            i = obj.nt + 1;
            obj.indices(i,:) = I;
            
            v1 = obj.vertices(:,I(1));
            v2 = obj.vertices(:,I(2));
            v3 = obj.vertices(:,I(3));
            
            %plane equation coefficients
            N = cross3(v2-v1,v3-v1);
            N = N/norm(N);
            obj.pec(i,:) = [N; -N'*v1];
            
            %inequality matrices
            if N(3) > 0
                obj.M12(i,:) = TriMeshSurf.getLineEqn(v1,v2);
                obj.M23(i,:) = TriMeshSurf.getLineEqn(v2,v3);
                obj.M31(i,:) = TriMeshSurf.getLineEqn(v3,v1);
            else
                %ignore these triangles during collision check
                obj.M12(i,:) = NaN;
                obj.M23(i,:) = NaN;
                obj.M31(i,:) = NaN;
            end

        end
        function addTriangles(obj, I)
            %I:     nx3, indices
            for i = 1:size(I,1)
                addTriangle(obj, I(i,:))
            end
        end
        function addQuad(obj, I)
            %I:     size 4, indices
            addTriangle(obj, I([1 2 3]))
            addTriangle(obj, I([1 3 4]))
        end
        

        function h = surfaceHeight(obj,pts)
            
            np = size(pts,2);
            h = NaN(1,np);

            loc = getLoc(obj,pts);
            in = loc > 0;
            
            if any(in)
                %z = -1/c*[a b d]*[x y 1]'
                pts_ = pts([1 2],in);
                pts_(3,:) = 1;
                h(in) = sum(obj.pec(loc(in),[1 2 4])' .* pts_,1);
                h(in) = -1./obj.pec(loc(in),3)' .* h(in);
            end

        end
        
        function dz = surfaceDz(obj,pts)
            np = size(pts,2);
            dz = NaN(1,np);
            
            loc = getLoc(obj,pts);
            in = loc > 0;
            if any(in)
                pts_ = pts(:,in);
                pts_(4,:) = 1; %append row of ones
                dz(in) = sum(obj.pec(loc(in),:)' .* pts_,1);
            end
            
        end

        function N = surfaceNormal(obj,pts)
            %TODO, add loc input?
            
            np = size(pts,2);
            N = NaN(3,np);
            
            loc = getLoc(obj,pts);
            in = loc > 0;
            
            N(:,in) = obj.pec(loc(in),1:3)';
            
        end
    end
    methods (Access = private)
        function loc = getLoc (obj, pts)
            %pts:   3 x np,
            %loc:   1 x np, index of triangle that each point is within, 0 if not within any triangle
            %TODO, can this be faster?
            
            np = size(pts,2);
            loc = zeros(1,np);
            for i = 1:np %loop over pts
                pt = pts(:,i);
                
                if pt(1) >= obj.xmin && pt(1) <= obj.xmax && ...
                   pt(2) >= obj.ymin && pt(2) <= obj.ymax
                    %in bounding box, now check if in any triangles
                    
                    pt(3) = 1;
                    in = obj.M12*pt <= 0 & obj.M23*pt <= 0 & obj.M31*pt <= 0; %nt x 1
                    loc_ = find(in,1,'first');
                    if ~isempty(loc_)
                        loc(i) = loc_;
                    end
                end
            end
        end
    end
    methods (Static)
        function lec = getLineEqn(v1,v2)
            %v1,v2:     3x1 vertices, only x y elements used
            %lec:       1x3 line equation coefficients for the line through the vertices [a b c]*[x y 1]' = 0
            
            a =  v2(2) - v1(2); %dy
            b =-(v2(1) - v1(1)); %-dx
            
            %[a b]' is a unit normal vector
            n = sqrt(a*a + b*b);
            a = a/n; b = b/n;

            c = -[a b]*v1(1:2);
            
            lec = [a b c]';
        end
        
    end
    
        
    
end
    
    