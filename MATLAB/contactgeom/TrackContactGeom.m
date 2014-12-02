classdef TrackContactGeom
    properties (GetAccess = 'public', SetAccess = 'public')
        dz          %1 x np
        HT_track    %4 x 4 x np, HT_contact_to_track
        HT_world    %4 x 4 x np, HT_contact_to_world
        dvc_dvt     %3*np x 7, Jacobian of contact point linear vel wrt [track frame spatial vel; sprocket angular rate]
        incontact   %1 x np, which points are in contact
    end
    
    properties (Dependent)
        np          %number of points
    end
    
    methods
        function out = get.np(obj)
            out = size(obj.HT_track,3);
        end
        
    end
end