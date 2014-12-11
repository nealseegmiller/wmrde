%controller function input/output

classdef ControllerIO
    properties (GetAccess = 'public', SetAccess = 'public')
        cmd         %na x 1, commanded
        interr      %na x 1, integrated error
        err         %na x 1, error: commanded - actual
        vis_act     %size nv logical, which elements of qvel are actuated. sum(vis_act) = na
    end

end