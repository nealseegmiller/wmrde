

classdef ForwardDynOutput < handle %handle class to avoid unnecessary copying in function calls
    %use copy() if necessary

    properties (GetAccess = 'public', SetAccess = 'public')
        vc0                 %3 x np, velocity of contact points wrt ground at time i
        vc                  %at time i+1
        modelf_contact      %3 x np, wheel/track-ground contact forces in contact coords
        modelf_act          %na x 1, actuator forces
        modelf_joint        %njc x 1, joint constraint forces
        %TODO, lambda?
        iter                %number of Newton's method iterations required for force-balance optimization
        cost                %final cost in force-balance optimization
    end
    

end