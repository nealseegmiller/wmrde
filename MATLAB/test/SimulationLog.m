%SimulationLog class


classdef SimulationLog < handle %handle class to avoid unnecessary copying in function calls
    %use copy() if necessary

    properties (GetAccess = 'public', SetAccess = 'private')

        time        %(n x 1)
        state       %(n x ns)
        qvel        %(n x nv) joint space vel
        contacts    %(1 x n) cell array of 
                    %1 x nw arrays of WheelContactGeom objects, or
                    %1 x nt arrays of TrackContactGeom objects
        u           %(n x na) actuated joint rate commands
        vc          %(3 x np x n) contact point velocities
        
        %for animation
        HT_parent   %(1 x n) cell array of (4 x 4 x nf) matrices
        
        %dynamic sim only
        qacc        %(n x nv) joint space acceleration
        fc          %(3 x np x n) contact forces
        fact        %(n x na) actuator forces
        iter        %(n x 1)
        cost        %(n x 1)

    end
    
    properties (Dependent)
        orient
        pos
        jd
        
        %contacts, TODO
        dz              %(n x np) contact height error
        contact_angles  %(n x np)

        
        %DEBUGGING
        vc_point            %permute such that n x 3 x np, easier to plot by point
        vc_dir              %permute such that n x np x 3, easier to plot by direction
        fc_point
        fc_dir
        
    end
    
    methods
        function obj = SimulationLog(n,nf,np,na,do_dyn)
            dlen = -1 + SIZEORIENT() + 3;
            ns = nf + dlen;
            nv = nf + 5;
            
            obj.time = NaN(n,1);
            obj.state = NaN(n,ns);
            obj.qvel = NaN(n,nv);
            obj.contacts = cell(1,n);
            obj.u = NaN(n,na);
            obj.vc = NaN(3,np,n);
            
            obj.HT_parent = cell(1,n);
            obj.HT_parent(:) = {zeros(4,4,nf)};
            
            if do_dyn %dynamic simulation
                obj.qacc = NaN(n,nv);
                obj.fc = NaN(3,np,n);
                obj.fact = NaN(n,na);
                obj.iter = NaN(n,1);
                obj.cost = NaN(n,1);
            end
        end
        function set_time(obj,i,in)
            obj.time(i) = in;
        end
        function set_state(obj,i,in)
            obj.state(i,:) = in;
        end
        function set_qvel(obj,i,in)
            obj.qvel(i,:) = in;
        end
        function set_contacts(obj,i,in)
            obj.contacts{i} = in;
        end
        function set_u(obj,i,in)
            obj.u(i,:) = in;
        end
        function set_vc(obj,i,in)
            obj.vc(:,:,i) = in;
        end
        function set_HT_parent(obj,i,in)
            obj.HT_parent{i}(:) = in;
        end
        function set_qacc(obj,i,in)
            obj.qacc(i,:) = in;
        end
        function set_fc(obj,i,in)
            obj.fc(:,:,i) = in;
        end
        function set_fact(obj,i,in)
            obj.fact(i,:) = in;
        end
        function set_iter(obj,i,in)
            obj.iter(i) = in;
        end
        function set_cost(obj,i,in)
            obj.cost(i) = in;
        end
        
        %DEPENDENT
        function out = get.orient(obj)
            isorient = get_stateIs(obj);
            out = obj.state(:,isorient);
        end
        function out = get.pos(obj)
            [~,ispos] = get_stateIs(obj);
            out = obj.state(:,ispos);
        end
        function out = get.jd(obj)
            [~,~,isjd] = get_stateIs(obj);
            out = obj.state(:,isjd);
        end
        
        function out = get.dz(obj)
            n = length(obj.contacts);
            np = size(obj.vc,2);
            out = NaN(n,np);
            for i = 1:n
                if ~isempty(obj.contacts{i})
                    out(i,:) = [obj.contacts{i}.dz];
                end
            end
        end
        function out = get.contact_angles(obj)
            n = length(obj.contacts);
            nw = length(obj.contacts{1});
            out = NaN(n,nw);
            for i = 1:n
                if ~isempty(obj.contacts{i})
                    out(i,:) = [obj.contacts{i}.angle];
                end
            end
        end

        
        function out = get.vc_point(obj)
            out = permute(obj.vc,[3 1 2]);
        end
        function out = get.vc_dir(obj)
            out = permute(obj.vc,[3 2 1]);
        end
        function out = get.fc_point(obj)
            out = permute(obj.fc,[3 1 2]);
        end
        function out = get.fc_dir(obj)
            out = permute(obj.fc,[3 2 1]);
        end
        

        %internal
        function [isorient,ispos,isjd] = get_stateIs(obj)
            ns = size(obj.state,2);
            [isorient,ispos,isjd] = stateIs(ns);
        end
        
        function [s,alpha] = get_calcSlip(obj,Rw,method)
            %Rw:        n x np, wheel radius * angular rate (get from .qvel)
            %method:    1,2
            
            vc_dir = obj.vc_dir;
            vx = vc_dir(:,:,1);
            vy = vc_dir(:,:,2);
            
            [s,alpha] = calcSlip(vx,vy,Rw,method);
        end
        
        function out = get_contact_angles_rot0(obj,rot)
            %rot:       n x nw, wheel rotation (get from .state)
            
            out = diffRad(obj.contact_angles,-rot);
        end
        
    end
    
end