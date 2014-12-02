classdef PlaneSurf < Surface
    
    properties (GetAccess = 'public', SetAccess = 'private')
        pec         %1x4, plane equation coefficients [a b c d]*[x y z 1]' = 0
    end

    methods
        function obj = PlaneSurf(in1,in2)
            
            obj.pec = zeros(1,4);
            
            %constructor
            if nargin > 1
                N = in1;    %normal vector, 3x1
                pt = in2;   %point, 3x1
                N = N/norm(N); %normalize
                obj.pec(:) = [N; -N'*pt];
            else
                obj.pec(:) = in1;
                N = obj.pec(1:3);
                N = N/norm(N); %normalize
                obj.pec(1:3) = N;

            end
            
        end

        %functions required by the Abstract superclass 'Surface'
        
        function h = surfaceHeight(obj,pts)
            %z = -1/c*[a b d]*[x y 1]'
            
            pts(3,:) = 1;
            h = -1/obj.pec(3) * obj.pec([1 2 4])*pts;
            
        end
        
        function dz = surfaceDz(obj,pts)
            %dz = [a b c d]*[x y z 1]'
            
            pts(4,:) = 1;
            dz = obj.pec*pts;
        end

        function N = surfaceNormal(obj,pts)
            n = size(pts,2);
            N = obj.pec(1:3)'*ones(1,n); %faster than repmat
        end
    end

end