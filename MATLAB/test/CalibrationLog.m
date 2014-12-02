%CalibrationLog class


classdef CalibrationLog < handle %handle class to avoid unnecessary copying in function calls
    %use copy() if necessary

    properties (GetAccess = 'public', SetAccess = 'private')
        %np:    size of parameter vector
        %nres:  size of entire residual (including stochastic)
        %nt:    number of residual types (e.g. calibrated, uncalibrated, no overlap)
        %ni:    number of iterations
        
        x0      %(np x 1) initial parameter values
        P0      %(np x np) initial parameter cov
        x       %(ni x np) parameter values after each iter.
        P       %(np x np x ni) parameter cov after each iter.
        
        I       %(ni x 2) interval indices [ind0 indf]
        res     %(ni x nres) matrix
        C       %(nres x nres x ni) matrix, residual covariances
        odomlog %(1 x ni) array of OdometryLog objects
        
        %DEBUGGING
        mdist   %(ni x 1) Mahalanobis distance
        md_tol  %scalar, Mahalanobis distance tolerance, if < 1 this is probability level such that md_tol = chi2inv(1-plevel,nres)
        
    end
    
    properties (Dependent)
        thresh_md
        thresh_plevel
        x0x
        P0P
        P0P_diag
    end
    
    
    methods
        function obj = CalibrationLog(np,nres,ni)
            obj.x0 = NaN(np,1);
            obj.P0 = NaN(np,np);
            obj.x = NaN(ni,np);
            obj.P = NaN(np,np,ni);
            
            obj.I = NaN(ni,2);
            obj.res = NaN(ni,nres);
            obj.C = NaN(nres,nres,ni);
            obj.odomlog = OdometryLog.empty(1,0);
            
            obj.mdist = NaN(ni,1);
        end
        
        function set_x0(obj,in)
            obj.x0(:) = in;
        end
        function set_P0(obj,in)
            obj.P0(:) = in;
        end
        function set_x(obj,i,in)
            obj.x(i,:) = in;
        end
        function set_P(obj,i,in)
            obj.P(:,:,i) = in;
        end
        function set_I(obj,i,in)
            obj.I(i,:) = in;
        end
        function set_res(obj,i,in)
            obj.res(i,:) = in;
        end
        function set_C(obj,i,in)
            obj.C(:,:,i) = in;
        end
        function set_odomlog(obj,i,in)
            obj.odomlog(i) = in;
        end
        function set_mdist(obj,i,in)
            obj.mdist(i) = in;
        end
        function set_md_tol(obj,in)
            obj.md_tol(1) = in;
        end
        
        %DEPENDENT
        function out = get.thresh_md(obj)
            if obj.md_tol < 1 %convert probability level to Mahalanobis distance
                nres = size(obj.res{1},2); %degrees of freedom
                out = chi2inv(1-obj.md_tol,nres);
            else
                out = obj.md_tol;
            end
        end
        function out = get.thresh_plevel(obj)
            if obj.md_tol < 1
                out = obj.md_tol;
            else
                nres = size(obj.res{1},2); %degrees of freedom
                out = 1-chi2cdf(obj.md_tol,nres);
            end
        end
        function out = get.x0x(obj)
            out = [obj.x0'; obj.x];
        end
        function out = get.P0P(obj)
            out = cat(3,obj.P0,obj.P);
        end
        function out = get.P0P_diag(obj)
            P0P = obj.P0P;
            n = size(P0P,3); %ni+1
            np = size(P0P,1);
            out = zeros(n,np);
            for i = 1:n
                out(i,:) = diag(P0P(:,:,i));
            end
        end
        
        
        function out = get_x_priorto(obj,ind0)
            %TODO, check this
            %get x calibrated only on data prior to ind0
            %assumes I(:,2) in increasing order
            
            i = find(obj.I(:,2) <= ind0,1,'last');
            if ~isempty(i)
                out = obj.x(i,:)';
            else
                out = obj.x0;
            end
        end
        
        function out = res_pos(obj,inres,inityaw0,inliers)
            %INPUTS
            %inres:         1 x ns logical, which elements of state are in residual
            %inityaw0:      if true, rotate such that initial yaw = 0, assumes Euler angles
            %inliers:       (optional) ni x 1 logical or []
            %OUTPUTS
            %out:           position residuals
            
            if nargin < 4
                inliers = [];
            end
            
            ni = size(obj.res,1);
            nr_ = sum(inres); %size of residual
            
            ns = get_ns(obj);
            [isorient,ispos] = get_stateIs(obj);
            
            out = NaN(ni,3);
            
            for i = 1:ni
                res_fullstate = zeros(ns,1);
                res_fullstate(inres) = obj.res(i,1:nr_);
                res_pos = res_fullstate(ispos);

                if inityaw0
                    euler = obj.odomlog(i).state(1,isorient);
                    R = Rotz(-euler(3));
                    res_pos = R*res_pos;
                end

                out(i,:) = res_pos;
            end
            
            if ~isempty(inliers)
                out = out(inliers,:);
            end
            
        end
        
        function out = C_pos(obj,inres,inityaw0,inliers)
            %out:           position residual covariances
            
            if nargin < 4
                inliers = [];
            end
            
            ni = size(obj.res,1);
            nr = sum(inres); %size of residual
            
            ns = get_ns(obj);
            [isorient,ispos] = get_stateIs(obj);
            
            out = NaN(3,3,ni);
            
            for i = 1:ni
                C_fullstate = zeros(ns,ns);
                C_fullstate(inres,inres) = obj.C(1:nr,1:nr,i);
                C_pos = C_fullstate(ispos,ispos);

                if inityaw0
                    euler = obj.odomlog(i).state(1,isorient);
                    R = Rotz(-euler(3));
                    C_pos = R*C_pos*R';
                end

                out(:,:,i) = C_pos;
            end
            
            if ~isempty(inliers)
                out = out(:,:,inliers);
            end
            
        end
        
        function out = res_yaw(obj,inres,inliers)
            %out:           yaw residuals
            
            if nargin < 3
                inliers = [];
            end
            
            ni = size(obj.res,1);
            nr = sum(inres); %size of residual
            
            ns = get_ns(obj);
            isorient = get_stateIs(obj);
            
            out = NaN(ni,1);
            
            for i = 1:ni
                res_fullstate = zeros(ns,1);
                res_fullstate(inres) = obj.res(i,1:nr);
                res_euler = res_fullstate(isorient);
                
                out(i) = res_euler(3);
            end
            
            if ~isempty(inliers)
                out = out(inliers);
            end
        end
        
        function out = C_yaw(obj,inres,inliers)
            %out:           yaw residual covariances
            
            if nargin < 3
                inliers = [];
            end
            
            ni = size(obj.res,1);
            nr = sum(inres); %size of residual
            
            ns = get_ns(obj);
            isorient = get_stateIs(obj);
            
            out = NaN(ni,1);
            
            for i = 1:ni
                C_fullstate = zeros(ns,ns);
                C_fullstate(inres,inres) = obj.C(1:nr,1:nr,i);
                C_euler = C_fullstate(isorient,isorient);
                
                out(i) = C_euler(3,3);
            end
            
            if ~isempty(inliers)
                out = out(inliers);
            end
        end
        
        function out = get_ns(obj)
            out = size(obj.odomlog(1).state,2);
        end
        function [out1,out2,out3] = get_stateIs(obj)
            ns = get_ns(obj);
            nf = size(obj.odomlog(1).qvel,2)-6+1;
            [out1,out2,out3] = stateIs(ns,nf);
        end
        
        function out = toStruct(obj)
            out = structNoDependent(obj);
%             out = struct(obj);
            
            %convert OdometryLog objects to structs
            out = rmfield(out,'odomlog');
            for i = 1:length(obj.odomlog)
                out.odomlog(i) = toStruct(obj.odomlog(i));
            end

        end
            
    end

    
end







