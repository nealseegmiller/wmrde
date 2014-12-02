classdef GridSurf < Surface

    properties (GetAccess = 'public', SetAccess = 'private')
        X
        Y
        Z %griddedInterpolant object, terrain height
        
        dx %X spacing
        dy %Y spacing
        
        %to compute the surface normal n = [nx ny nz]';
        DZDX %griddedInterpolant object, gradient wrt x
        DZDY %griddedInterpolant object, gradient wrt y
        
    end
    methods
        function obj = GridSurf(xmin,xmax,nx,ymin,ymax,ny,Z)
            %X,Y matrices from inputs
            %griddedInterpolant requires ndgrid
            [X_,Y_] = ndgrid(linspace(xmin,xmax,nx),linspace(ymin,ymax,ny));
            
            if any(size(X_) ~= size(Z)) || any(size(Y_) ~= size(Z))
                disp('MeshSurf: size of X or Y does not match size of Z')
                return
            end
            
            obj.X=X_;
            obj.Y=Y_;
            obj.Z=griddedInterpolant(X_,Y_,Z,'linear');

            
            obj.dx = (xmax-xmin)/(nx-1);
            obj.dy = (ymax-ymin)/(ny-1);
            
            %help gradient, 
            %1st output is differences in horizontal direction
            %2nd output is differences in vertical direction
            [DZDY_, DZDX_] = gradient(Z); 
            DZDX_ = DZDX_/obj.dx; %dz/dx
            DZDY_ = DZDY_/obj.dy; %dz/dy
            
            obj.DZDX=griddedInterpolant(X_,Y_,DZDX_,'linear');
            obj.DZDY=griddedInterpolant(X_,Y_,DZDY_,'linear');
            
            
        end
        
        
        %define the functions required by the Abstract superclass 'Surf'
        function h = surfaceHeight(obj,pts)
            h = obj.Z(pts(1,:),pts(2,:));
        end
        
        function dz = surfaceDz(obj,pts)
            z_surf = surfaceHeight(obj,pts(1:2,:));
            
            dh = (pts(3,:)-z_surf);
            
            %dz = dot product of delta height (along world z axis) and surface normal
            N = surfaceNormal(obj,pts(1:2,:)); %expensive!
            dz = dh.*N(3,:);
            
        end
        
        function N = surfaceNormal(obj,pts)
            n = size(pts,2);
            
            %interpolate normal is NOT equivalent to finite diff.
            dzdx = obj.DZDX(pts(1,:),pts(2,:));
            dzdy = obj.DZDY(pts(1,:),pts(2,:));
            
            N = [-dzdx; -dzdy; ones(1,n)];
            
            %normalize
            normN = sqrt(sum(N.^2,1));
            N = N./[normN; normN; normN];
            
        end
    end
    
end