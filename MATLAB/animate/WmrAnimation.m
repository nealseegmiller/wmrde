classdef WmrAnimation < handle %handle class to avoid unnecessary copying in function calls
    %use copy() if necessary

    properties (GetAccess = 'public', SetAccess = 'public')
        ax_center %3x1 (in body coordinates)
        ax_limits %1x6
        
        h_quiver %TODO, for track contact frames
    end
    
    properties (GetAccess = 'public', SetAccess = 'private')
        parent_ind %parent indices for hgtransforms
        h_fig
        h_axis
        h_light
        h_hgt %hgtransform, one for each frame
        h_line = cell(1,10); %cell array
        
    end

    methods
        function obj = WmrAnimation(FigureHandle,AxisHandle)
            
            if nargin > 0
                %use existing figure and axis
                obj.h_fig = FigureHandle;
                obj.h_axis = AxisHandle;
                
            else
                set(figure,'name','WMR animation');
                obj.h_fig = gcf;
                obj.h_axis = gca;
            
                hold on

                set(gcf,'Renderer','OpenGL')
%                 set(gcf,'Renderer','zbuffer') %slower

                xlabel('x')
                ylabel('y')
                zlabel('z')
                axis equal

                obj.h_light = light('Position',[1 0 1],'Parent',obj.h_axis);
                
                makeLegible(14)
            end

        end
        
        function setHgtParent(obj,index,parent_index)
            if parent_index > 0
                set(obj.h_hgt(index),'Parent',obj.h_hgt(parent_index))
            else %<0 or NaN, body frame
                set(obj.h_hgt(index),'Parent',obj.h_axis) 
            end
            obj.parent_ind(index) = parent_index;
        end
        
        function addHgt(obj,parent_index)   
            for i = 1:length(parent_index)
                obj.h_hgt(end+1) = hgtransform();
                setHgtParent(obj,length(obj.h_hgt),parent_index(i));
            end
        end
        
        function loadHgt(obj,filename)
            obj.h_hgt = hgload(filename);
            load([filename '_parentind.mat'])
            for i = 1:length(obj.h_hgt)
                setHgtParent(obj,i,ParentInd(i));
                set(obj.h_hgt(i),'Matrix',eye(4)); %DEBUGGING
            end
        end
        
        function saveHgt(obj,filename)
            %hgsave requires that handles belong to unrelated objects
            %Parent everything to axis, then put it back
            ParentInd = obj.parent_ind; %backup
            
            for i = 1:length(obj.h_hgt)
                setHgtParent(obj,i,NaN);
            end
            hgsave(obj.h_hgt,filename)
            save([filename '_parentind.mat'],'ParentInd')
            for i = 1:length(obj.h_hgt)
                setHgtParent(obj,i,ParentInd(i));
            end
        end
            
        
        function updateHgt(obj,HT_parent,inds)
            %update hgtransform matrices
            %HT_parent:	4x4xn homogeneous transforms to parent coords
            %inds:      (optional) 1xn vector of hgtransform indices
            
            if nargin < 3
                inds = 1:length(obj.h_hgt);
            end
            
            for i = 1:length(inds)
                set(obj.h_hgt(inds(i)),'Matrix',HT_parent(:,:,i))
            end
        end
        
        function centerAxes(obj,HT_body_to_world)
            c = applyHT(HT_body_to_world,obj.ax_center);
            C = [c(1) c(1) c(2) c(2) c(3) c(3)];
            
            axis(C + obj.ax_limits);

        end
        
        function h=addQuiverZ(obj,inds,len)
            %add quiver3 on z axis of hgtransform
            %inds:	1xn vector of indices
            %len:   scalar, quiver length
            %OUTPUT
            %h:     1xn vector of handles
            
            xyz = [0 0 0]';
            uvw = [0 0 1]'*len;
            h = zeros(1,length(inds));
            for i = 1:length(inds)
                h(i) = drawQuivers(xyz,uvw);
                set(h(i),'Parent',obj.h_hgt(inds(i)));
            end
        end
        
        function h=addQuiverXYZ(obj,inds,len)
            %add quiver3 on x,y,z axes of hgtransform
            %inds:	1xn vector of indices
            %len:   scalar, quiver length
            %OUTPUT
            %h:     1xn vector of handles
            
            xyz = zeros(3);
            uvw = eye(3)*len;
            h = zeros(1,length(inds));
            for i = 1:length(inds)
                h(i) = drawQuivers(xyz,uvw);
                set(h(i),'Parent',obj.h_hgt(inds(i)));
            end
        end
        
        function h=addText(obj,inds,str)
            %add text to hgtransform
            %inds:  1xn vector of indices
            %str:	string or size n cell array of strings
            
            if ~iscell(str), str = {str}; end
            
            h = zeros(1,length(inds));
            for i = 1:length(inds)
                h(i) = str(0,0,0,str{i});
                set(h(i),'Parent',obj.h_hgt(inds(i)));
            end
        end
        
        
        function updateLine(obj,setno,pts)
            %append a point to each line in set, create line if necessary
            %setno:	scalar, line set number
            %pts:   3xn points, n equals number of lines in set


            if isempty(obj.h_line{setno})
                %make new
                for i = 1:size(pts,2)
                    obj.h_line{setno}(i) = plot3(pts(1,i),pts(2,i),pts(3,i));
                    set(obj.h_line{setno}(i),'Clipping','off')
                end
            else
                %append
                for i = 1:length(obj.h_line{setno})
                    X = get(obj.h_line{setno}(i),'XData');
                    Y = get(obj.h_line{setno}(i),'YData');
                    Z = get(obj.h_line{setno}(i),'ZData');

                    X(end+1) = pts(1,i); %#ok<AGROW>
                    Y(end+1) = pts(2,i); %#ok<AGROW>
                    Z(end+1) = pts(3,i); %#ok<AGROW>

                    set(obj.h_line{setno}(i),'XData',X);
                    set(obj.h_line{setno}(i),'YData',Y);
                    set(obj.h_line{setno}(i),'ZData',Z);
                end
            end
            
        end
        
        function  deleteLine(obj,setno)
            delete(obj.h_line{setno})
            obj.h_line{setno} = [];
        end
        
        
        function updateAnimation(obj, HT_to_parent, contacts)
            %modify this function as needed
            
            init = true; %flag, true if initialized
            if isempty(obj.h_line{1})
                init=false;
            end

            nf = size(HT_to_parent,3);
            updateHgt(obj, HT_to_parent, 1:nf)
            
            if isa(contacts,'WheelContactGeom')
                %contact frames
                
                nw = length(contacts);
                I = (nf+1):(nf+nw);
                updateHgt(obj, cat(3,contacts.HT_wheel), I)
                
                %visible only if in contact
                for wno = 1:nw
                    vis = 'off';
                    if contacts(wno).incontact
                        vis = 'on';
                    end
                    set(obj.h_hgt(I(wno)),'Visible',vis);
                end
                
                cp = contactPointsOnSurface(contacts);
                updateLine(obj, 2, cp); %contact points
                
                
                
            elseif isa(contacts,'TrackContactGeom')
                %track frames
                
                nt = length(contacts);
                
                for tno = 1:nt
                    np = size(contacts(tno).HT_track,3);
                    cp = zeros(3,np);
                    cp(:) = contacts(tno).HT_track(1:3,4,:);
                    
                    N = zeros(3,np);
                    N(:) = contacts(tno).HT_track(1:3,3,:);
                    
                    
                    incontact = [contacts(tno).incontact];
                    d = .05; %quiver length
                    
%                     %refresh error if number of quivers changes
%                     drawQuivers(cp(:,incontact),N(:,incontact)*d,obj.h_quiver(tno)); 
                    
                    uvw = N .*(d*ones(3,1)*incontact);
                    drawQuivers(cp,uvw,obj.h_quiver(tno)); 
                    
                end
                
            end
            
            
            updateLine(obj, 1, HT_to_parent(1:3,4,1)); %body
            
            
            if ~init
                set(obj.h_line{1},'LineWidth',1.5,'Color','b')
                set(obj.h_line{2},'LineWidth',1.5,'Color','r')
            end
            
            centerAxes(obj, HT_to_parent(:,:,1));
            drawnow
            
        end
        
        function reset(obj)
            axes(obj.h_axis);
            
            %delete all lines
            for i = 1:length(obj.h_line)
                if ~isempty(obj.h_line{i})
                    deleteLine(obj,i);
                end
            end
            
            %set hgtransform matrices to identity
            HT_to_parent = repmat(eye(4),[1 1 length(obj.h_hgt)]);
            updateHgt(obj, HT_to_parent);
            
            centerAxes(obj, HT_to_parent(:,:,1));
            
            drawnow
        end
        
    end
    
%     methods (Static)
% 
%     end
    
end