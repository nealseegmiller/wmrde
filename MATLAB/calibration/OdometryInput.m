classdef OdometryInput < handle
    
    properties (GetAccess = 'public', SetAccess = 'public')
        %nv: length of joint space velocity
        %nm: number of measurements
        %ns: length of state vector
        t       %n x 1, time
        u       %n x na, controller inputs, NaN for free elements
        z       %n x nm, measurements
                %may also be n x ns, NaN for not-measured elements of state
        H       %nm x ns, measurement Jacobian
        R       %nm x nm, measurement noise
        Q       %nv x nv, spectral density of joint space velocity (without added noise to body vel.)
        %TODO, let z, H, R, Q vary with time. requires cell arrays if nm varies
        
        vis_act %size nv logical, which elements of qvel are actuated. sum(vis_act) = na
    end
end