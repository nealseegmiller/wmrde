%OdometryLog class


classdef OdometryLog < handle %handle class to avoid unnecessary copying in function calls
    %use copy() if necessary

    properties (GetAccess = 'public', SetAccess = 'private')

        time    %(n x 1)
        state   %(n x ns)
        qvel    %(n x nv)
        %(optional)
        P       %(ns x ns x n) %state covariance
        Q       %(nv x nv x n) %qvel spectral density
        
        HT_sensor_to_body = eye(4) %(4x4) homogeneous transform from sensor to body coords

        %DEBUGGING
        vc                  %(3 x nw x n) contact point velocities
        contacts            %1 x n cell array of 1 x nw WheelContactGeom objects
        iter                %(n x 1)
    end
    
    properties (Dependent)        
        pos
        pos_sensor
        vel_sensor
        P_diag %TODO
        
        %DEUBGGING
        vc_wheel            %permute such that n x 3 x nw, easier to plot by wheel
        vc_dir              %permute such that n x nw x 3, easier to plot by direction
    end
    
    methods
        function obj = OdometryLog(n,ns,nv,nw,do_stoch)
            obj.time = NaN(n,1);
            obj.state = NaN(n,ns);
            obj.qvel = NaN(n,ns);

            if do_stoch
                obj.P = NaN(ns,ns,n);
                obj.Q = NaN(nv,nv,n);
            end

            %DEBUGGING
            obj.vc = NaN(3,nw,n);
            obj.contacts = cell(1,n);
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
        function set_P(obj,i,in)
            obj.P(:,:,i) = in;
        end
        function set_Q(obj,i,in)
            obj.Q(:,:,i) = in;
        end
        function set_HT_sensor_to_body(obj,in)
            obj.HT_sensor_to_body(:) = in;
        end
        %DEBUGGING
        function set_vc(obj,i,in)
            obj.vc(:,:,i) = in;
        end
        function set_contacts(obj,i,in)
            obj.contacts{i} = in;
        end
        function set_iter(obj,i,in)
            obj.iter(i) = in;
        end
        
        %sensor
        function out = get.pos(obj)
            ns = size(obj.state,2);
            nf = size(obj.qvel,2)-5;
            [~,ispos] = stateIs(ns,nf);
            out = obj.state(:,ispos);
        end
        
        function out = get.pos_sensor(obj)
            %output position of sensor in world coordinates
            n = size(obj.state,1);
            ns = size(obj.state,2);
            nf = size(obj.qvel,2)-5;
            [isorient,ispos] = stateIs(ns,nf);
            tsensor = obj.HT_sensor_to_body(1:3,4);
            
            out = zeros(n,3);
            for i = 1:n
                HT_body_to_world = poseToHT(obj.state(i,isorient)',obj.state(i,ispos)');
                out(i,:) = applyHT(HT_body_to_world,tsensor)';
            end
        end
        function out = get.vel_sensor(obj)
            %output spatial velocity of sensor in body coordinates
            Xm = invHTToPlucker(obj.HT_sensor_to_body);
            out = (Xm*obj.qvel(:,1:6)')';
        end
        function out = get.P_diag(obj)
            %diagonals of P matrices
            n = size(obj.P,3);
            ns = size(obj.P,1);
            out = zeros(n,ns);
            for i = 1:n
                out(i,:) = diag(obj.P(:,:,i));
            end
        end
        function out = get.vc_wheel(obj)
            out = permute(obj.vc,[3 1 2]);
        end
        function out = get.vc_dir(obj)
            out = permute(obj.vc,[3 2 1]);
        end
        
        function [h,h_end,he_line,he_surf,he_cent] = plotPath(obj,do_sensor,I,dim,marg,h_axis)

            if nargin < 6
                h_axis = []; %existing axis handle
                if nargin < 5
                    marg = 3; %for covariance ellipse in 3D plot only, {1,2,3,[]}
                    if nargin < 4
                        dim = 2; %2D or 3D plot
                        if nargin < 3
                            I = []; %indices at which to draw covariance ellipse
                            if nargin < 2
                                do_sensor = false; %plot sensor position
                            end
                        end
                    end
                end
            end
            
            if isempty(h_axis)
                set(figure,'name','path')
                hold on
                xlabel('x')
                ylabel('y')
                zlabel('z')
            else
                axes(h_axis)
            end
            
            if do_sensor
                pos = obj.pos_sensor;
            else
                pos = obj.pos;
            end
            
            h_end = zeros(2,1);
            if dim == 2
                h = plot(pos(:,1),pos(:,2));
                h_end(1) = plot(pos(1,1),pos(1,2),'>'); %start
                h_end(2) = plot(pos(end,1),pos(end,2),'s'); %stop
                marg = 3;
            elseif dim == 3
                h = plot3(pos(:,1),pos(:,2),pos(:,3));
                h_end(1) = plot3(pos(1,1),pos(1,2),pos(1,3),'>'); %start
                h_end(2) = plot3(pos(end,1),pos(end,2),pos(end,3),'s'); %stop
            else
                disp([mfilename ': invalid dimension'])
            end

            he_line = [];
            he_surf = [];
            he_cent = [];
            if ~isempty(I) && ~isempty(obj.P)
                nvertices = 36*2;
                conf = .95;
                
                for i = I(:)'
                    ns = size(obj.state,2);
                    nf = size(obj.qvel,2)-5;
                    [~,ispos] = stateIs(ns,nf);
                    
                    [hl,hs] = drawEllipse(obj.P(ispos,ispos,i),pos(i,:)',nvertices,conf,marg);
                    he_line = [he_line; hl];
                    he_surf = [he_surf; hs];
                    he_cent = [he_cent; plot3(pos(i,1),pos(i,2),pos(i,3),'.')]; %center
                end

            end
            
            if isempty(h_axis)
                axis equal
                makeLegible(14)
            end
            
        end
        
        function out = toStruct(obj)
            out = structNoDependent(obj);
%             out = struct(obj);
        end

    end
    
end