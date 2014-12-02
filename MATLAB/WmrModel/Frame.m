%member of WmrModel class

classdef Frame %NOT a handle class!
    
    properties (GetAccess = 'public', SetAccess = 'public')
        name            %string
        dof_type        %1-6
        parent_ind      %index of parent frame
        HT_parent_jd0   %4x4 Homogeneous transform to parent coords if joint displacement == 0
        
        %logicals:
        iswheel
        issprocket      %for track
        % iswheel and issprocket cannot both be true
        isactuated
        isfixed         %if true, keep fixed when initializing terrain contact, == isactuated by default
        
        %for wheel & sprocket frames only
        rad             %radius of wheel or sprocket
        
%         dz_target       %TODO
%         tc_z
%         erp_z
%         cfm_z
%         fds_x
%         fds_y

        %for sprocket frames only
        rad2            %radius of 2nd sprocket
        L               %distance between sprockets
        
%         min_npic        %minimum number of points in contact, for kinematic sim.

    end
    
    properties (GetAccess = 'public', SetAccess = 'private')
        %mass properties
        %TODO
        %track mass is affixed to sprocket's parent frame
        %sprocket frame has no mass, but rotational inertia of entire track
        
        mass
        cm              %3x1, center of mass
        I               %3x3, moment of inertia
        
        %keep consistent
        %could be Dependent properties, but precompute for speed
        Is              %spatial inertia
    end
    
    methods
        function obj = Frame()

        end
        
        function obj = setMass(obj, scalar_mass, center_of_mass, moment_of_inertia)
            %[] to leave unchanged

            if ~isempty(scalar_mass)
                obj.mass = scalar_mass; 
            end
            if ~isempty(center_of_mass)
                obj.cm = center_of_mass;
            end
            if ~isempty(moment_of_inertia)
                obj.I = moment_of_inertia;
            end

            %update spatial inertia
            obj.Is = toSpatialInertia(obj.mass,obj.cm,obj.I);
        end
    end
    
    methods (Static)
        
        function dof_int = convertDofTypeString(dof_string)
            %convert string to integer
            %that specifies degree of freedom type
            switch dof_string
                case 'RX', dof_int = 1; %revolute
                case 'RY', dof_int = 2;
                case 'RZ', dof_int = 3;
                case 'PX', dof_int = 4; %prismatic
                case 'PY', dof_int = 5;
                case 'PZ', dof_int = 6;
                otherwise, dof_int = NaN;
            end
        end

    end
end