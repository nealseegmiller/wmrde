%additional outputs for odeKin, odeDyn functions

classdef OdeOutput < handle
    
    properties (GetAccess = 'public', SetAccess = 'public')
        %both
        HT_parent           %4x4x nf, homogeneous transform to parent coords
        HT_world            %4x4x nf , to world coords
        contacts            %1 x nw array of WheelContactGeom objects, or
                            %1 x nt array of TrackContactGeom objects
        u                   %na x 1, controller function output
        vc                  %3 x np, contact point velocities in contact coords
        motiontime          %scalar, computation time for motion prediction
        
        %odeDyn only
        fc                  %3 x np, wheel/track-ground contact forces in contact coords
        fact                %na x 1, actuator forces
        iter
        cost

    end
    
end