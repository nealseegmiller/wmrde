%WmrModel class


classdef WmrModel < handle %handle class to avoid unnecessary copying in function calls
    %use copy() if necessary
    properties (GetAccess = 'public', SetAccess = 'private')

        nf              %number of frames
        frames=Frame.empty(1,0) %1xnf array of Frame objects
        
        %could be Dependent, but precomputed for speed
        nw
        wheelframeinds
        nt              %number of tracks = number of sprocket frames
        sprocketframeinds
        na
        actframeinds
        
        Is              %spatial inertia of entire WMR
    end
    
    properties (GetAccess = 'public', SetAccess = 'public')
        
        
        %TODO, make these private?
        description     %string
        
        %FOR KINEMATIC MODEL
        min_npic    %input to whichPointsInContact
        dz_target   %1 x nw, target Delta z (contact height error)
        %time constants for Baumgarte's stabilization method
        tc_z        %1 x nw, contact constraint, z dir
        tc_j        %1 x njc, holonomic joint constraints
        
        %FOR DYNAMIC MODEL
        grav = 9.81;        %scalar acceleration of gravity
        use_constraints = 1; %else use stiff ODE
                
        
        %FUNCTION HANDLES AND PARAMETERS
        controller_fh;
        %[u,cmd] = controller(obj, time, state)
        %INPUT
        %obj:       WmrModel object
        %time:      scalar
        %state:     ns x 1
        %OUTPUT
        %u:         na x 1, joint rates for all actuated joints
        %cmd:       nv x 1, commanded qvel
        
        bsm_fh; %body slip model function handle
        bsm_p; %body slip model parameters
        %s = feval(bsm_fh,bsm_p,vbody,g,Is)
        %INPUT
        %vbody: 6x1, input spatial velocity of body frame wrt world in body coords
        %g:     3x1, grav. acceleration in body coords
        %Is:    6x6, spatial inertia of WMR in body coords
        %OUTPUT
        %s:     6x1, spatial slip velocity of body frame
        
        wgc_fh; %wheel-ground contact model function handle
        wgc_p; %wheel-ground contact model parameters
        %[f,J] = feval(wgc_fh,wgc_p,Iw,vc,Rw,dz)
        %INPUT
        %Iw:    1 x nw, wheel indices
        %vc:    3 x nw, velocity of contact pt wrt ground (in contact coords), [vx; vy; vz]
        %Rw:    1 x nw, wheel radius * wheel angular velocity
        %dz:    1 x nw, contact height error
        %OUTPUT
        %f:     3 x nw, contact forces exerted on wheel by the ground (in contact coords), [fx;fy;fz] 
        %J:     (optional) 3x5 x nw, Jacobian
        %       J(:,:,i) is the Jacobian of [fx fy fz]' with respect to [vx vy vz Rw dz]' for wheel i
        
        act_fh; %actuator model function handle
        act_p; %actuator parameters
        %[tau,err,dtaudu] = feval(act_fh,act_p,ucmd,u,interr)
        
        hjc_fh; %function handle for holonomic joint constraints
        %[c, Jc, f, df_djd, df_djr] = feval(hjc_fh, mdl, jd, jr)
        %INPUT
        %mdl:       WmrModel object
        %jd:        nj x 1, joint displacements
        %jr:        nj x 1, joint rates (d/dt jd)
        %OUTPUT
        %c:         nc x 1, constraint function error
        %Jc:        nc x nj, Jacobian of c wrt joint displacements
        %f:         nc x 1, force
        %df_djd:    nc x nj, Jacobian of force wrt joint displacements
        %df_djr:    nc x nj, Jacobian of force wrt joint rates
        
        cov_p %parameters of body frame spatial velocity covariance (Q), for stochastic model
        
        
        %the following parameters are only used if no wheel-ground contact model function is provided
        
        %1 x nw, for wheel-ground contact constraints
        erp_z       %error reduction parameter, z dir
        cfm_z       %constraint force mixing, z dir
        fds_x       %force dependent slip, x dir
        fds_y       %force dependent slip, y dir
        %1 x njc, for holonomic joint constraints
        erp_j
        cfm_j
        
    end



    methods
        function obj = WmrModel(fb)
            %constructor
            obj.nf = 0;
            for i = 1:fb
                obj.frames(i) = Frame();
            end
            
            %Dependent
            obj.nw = 0;
            obj.na = 0;

        end
        
        function initFrame(obj,fi)
            %initialize frame
            
            obj.frames(fi).name = [];
            obj.frames(fi).dof_type = 0;
            obj.frames(fi).parent_ind = 0;
            obj.frames(fi).HT_parent_jd0 = eye(4);

            obj.frames(fi).iswheel = false;
            obj.frames(fi).issprocket = false;
            obj.frames(fi).isactuated = false;
            obj.frames(fi).isfixed = false;

            obj.frames(fi) = setMass(obj.frames(fi),0,zeros(3,1),zeros(3)); 
        end
        
        %FUNCTIONS TO ADD FRAMES TO TREE
        function addBodyFrame(obj,name)
            %ERROR CHECKING
            if ~ischar(name)
                disp([mfilename ': frame names must be strings'])
                return
            end
            
            %set
            i = 1;
            obj.nf = i;
            initFrame(obj,i);
            obj.frames(i).name = name;
            
        end
        
        
        function addJointFrame(obj,name,parent_name,dof_string,isactuated,HT_parent_jointdisp0)
            
            %ERROR CHECKING
            if obj.nf >= length(obj.frames)
                disp([mfilename ': Buffer exceeded, no more frames may be added'])
                return
            end
            
            if ~ischar(name) || ~ischar(parent_name)
                disp([mfilename ': frame names must be strings'])
                return
            end
            
            names = {obj.frames.name};
            if any(strcmp(name,names))
                disp([mfilename ': ''name'' invalid; already in use']);
                return
            end
            
            pfi = namesToInds(obj,parent_name); %parent frame index
            
            if pfi==0
                disp([mfilename ': ''parent_name'' invalid, must match existing frame'])
                return
            end
            
            if obj.frames(pfi).iswheel
                disp([mfilename ': ''parent_name'' invalid, parent cannot be wheel frame'])
                return
            end
            
            dof_int = Frame.convertDofTypeString(dof_string);
            if dof_int == 0
                disp([mfilename ': ''dof_string'' input invalid'])
                return
            end
            
            if ~isValidHT(HT_parent_jointdisp0)
                disp([mfilename ': ''HT_parent_jointdisp0'' input is invalid'])
                return
            end

            i = obj.nf+1;
            obj.nf = i;
            initFrame(obj,i);
            obj.frames(i).name = name;
            obj.frames(i).parent_ind = pfi;
            obj.frames(i).dof_type = dof_int;
            obj.frames(i).HT_parent_jd0 = HT_parent_jointdisp0;
            obj.frames(i).iswheel = false;
            obj.frames(i).issprocket = false;
            obj.frames(i).isactuated = isactuated;
            obj.frames(i).isfixed = isactuated;
            
            isactuated = [obj.frames.isactuated];
            obj.na = sum(isactuated);
            obj.actframeinds = find(isactuated);

        end
        
        function addWheelFrame(obj, name, parent_name, isactuated, HT_parent_jointdisp0, radius)            
            
            %by convention, wheels rotate about their y axis
            addJointFrame(obj, name, parent_name, 'RY', isactuated, HT_parent_jointdisp0);
            i = obj.nf;
            obj.frames(i).iswheel = true;
            obj.frames(i).rad = radius;
            
            iswheel = [obj.frames.iswheel];
            obj.nw = sum(iswheel);
            obj.wheelframeinds = find(iswheel);
                        
        end
        
        function addSprocketFrame(obj,name,parent_name,isactuated, HT_parent_jointdisp0, radius, radius2, len)
            %by convention, sprockets rotate about their y axis
            addJointFrame(obj, name, parent_name, 'RY', isactuated, HT_parent_jointdisp0);
            i = obj.nf;
            obj.frames(i).issprocket = true;
            obj.frames(i).rad = radius;
            obj.frames(i).rad2 = radius2;
            obj.frames(i).L = len;
            
            issprocket = [obj.frames.issprocket];
            obj.nt = sum(issprocket);
            obj.sprocketframeinds = find(issprocket);
        end

        
        %PROPERTY ACCESS METHODS
        
        function set_HT_parent_jd0(obj, fi, ri, ci, HT_parent_jointdisp0_)
            %INPUT
            %fi:    frame index
            %ri:    row indices, [] for all
            %ci:    col indices, [] for all
            
            if isnan(obj.frames(fi).dof_type)
                disp([mfilename ': cannot set HT_parent_jd0 for this frame'])
            end
            if isempty(ri)
                ri=1:4; 
            end
            if isempty(ci)
                ci=1:4; 
            end
            obj.frames(fi).HT_parent_jd0(ri,ci) = HT_parent_jointdisp0_;
        end
        
        function set_rad(obj, fi, Radius)
            if ~obj.frames(fi).iswheel
                disp([mfilename ': cannot set radius, not a wheel frame'])
            end
            obj.frames(fi).rad = Radius;
        end
        
        function setFrameMass(obj, fi, mass, center_of_mass, moment_of_inertia)
            obj.frames(fi) = setMass(obj.frames(fi), mass, center_of_mass, moment_of_inertia);
        end
        
        function set_isfixed(obj, fi, val)
            obj.frames(fi).isfixed = val;
        end
        
        function out = namesToInds(obj,names_in)
            %returns the indices of frames by name
            %names_in:  a string or size n cell array of strings
            %out:       1 x n vector of indices, 0 if no match
            
            if ~iscell(names_in)
                names_in = {names_in};
            end
            
            frame_names = {obj.frames.name};
            out = zeros(1,length(names_in));
            for i = 1:length(names_in)
                %return the first match (should be the only match)
                fi = find(strcmp(frame_names,names_in{i}),1);
                if isempty(fi)
                    fi=0;
                end
                out(i) = fi;
            end
            
        end
        
        function update_Is(obj,HT_parent)
            if nargin < 2
                HT_parent = cat(3,obj.frames.HT_parent_jd0);
            end
            %total spatial inertia
            Xup = invHTsToPluckers(HT_parent);
            tmp = subtreeInertias(obj,Xup);
            obj.Is = tmp(:,:,1); 
        end
        
    end
    
    
end




