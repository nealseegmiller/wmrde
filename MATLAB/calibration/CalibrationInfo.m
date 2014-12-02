
classdef CalibrationInfo %Calibration Information

    properties (GetAccess = 'public', SetAccess = 'public')
        %np:    number of parameters
        %ns:    number of states
        param_names     %(1 x np cell array of strings)
        isstoch         %(1 x np logical) is stochastic parameter
        isfixed         %(1 x np logical) is fixed parameter
        inresidual      %(1 x ns logical) which elements of state are in the residual
    end
    
end