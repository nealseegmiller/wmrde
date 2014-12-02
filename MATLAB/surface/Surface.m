classdef Surface < handle %make a handle class to avoid unnecessary copying of mesh data
    %Matlab requires that:
    %If a class defines super-classes, all or none must be handle classes.
    %i.e. can't do 'classdef GridSurf < Surface & handle' if Surf is a value class
    %this means *all* subclasses of Surface must be handle classes
    
    methods (Abstract)
        %these methods must be implemented by subclasses
        %the inputs should be the same in all subclasses but Matlab doesn't enforce this
        
        h = surfaceHeight(obj,pts)
        %compute surface height (world z axis)
        %INPUT
        %pts (2xn matrix) contains x,y point locations in world coords
        %OUTPUT
        %h (1xn vector) height of terrain surface at the specified points
        
        dz = surfaceDz(obj,pts)
        %compute the distance of points from the surface (measured along the surface normal)
        %sign convention: negative for points below surface
        %INPUT
        %pts (3xn matrix) contains x,y,z point locations in world coords
        %OUTPUT
        %dz (1xn vector)
            
        N = surfaceNormal(obj,pts)
        %compute the surface normals
        %INPUT
        %pts (2xn matrix) contains x,y, point locations in world coords
        %OUTPUT
        %N (3xn matrix) each column is a normal vector
        
        %when 2xn pts matrix is required, 3xn pts matrix can be accepted but 3rd row is ignored
        
    end
   
end