classdef WheelContactGeom
    properties (GetAccess = 'public', SetAccess = 'public')
        dz          %contact height error
        angle       %contact angle
        HT_wheel    %4x4, HT_contact_to_wheel
        HT_world    %4x4, HT_contact_to_world
        incontact
        np = 1;
    end
end