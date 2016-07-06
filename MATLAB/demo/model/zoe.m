function [mdl, state, qvel, anim] = zoe()

opts.fix_front_axle_roll = 1;

%Dimensions
L = 1.91;
B = 1.64;
D = .119;
%wheels
rad = .325; %wheel radius

% B = 1.70; %DEBUGGING
nw=4; %number of wheels
%Masses
TotalMass=198; %kg
Mb = .9*TotalMass;
Mw = .1*TotalMass/nw;

Hb = .5; %height of body
Wb = 1; %width of body
Ww = .07; %width of wheel


mdl = WmrModel(8+~opts.fix_front_axle_roll);
addBodyFrame(mdl,'body');

%BUILD WMR MODEL OBJECT
if opts.fix_front_axle_roll
    addJointFrame(mdl,'s1','body','RZ',0,poseToHT([0 0 0],[L/2, 0, 0]));
else
    addJointFrame(mdl,'r1','body','RX',0,poseToHT([0 0 0],[L/2, 0, 0]));
    addJointFrame(mdl,'s1','r1','RZ',0,poseToHT([0 0 0],[0, 0, 0]));
end

addWheelFrame(mdl,'FL','s1',1,poseToHT([0 0 0],[0, B/2,-D]),rad);
addWheelFrame(mdl,'FR','s1',1,poseToHT([0 0 0],[0,-B/2,-D]),rad);

%rear axle
addJointFrame(mdl,'r2','body','RX',0,poseToHT([0 0 0],[-L/2, 0, 0]));
addJointFrame(mdl,'s2','r2','RZ',0,poseToHT([0 0 0],[0, 0, 0]));

addWheelFrame(mdl,'BL','s2',1,poseToHT([0 0 0],[0, B/2,-D]),rad);
addWheelFrame(mdl,'BR','s2',1,poseToHT([0 0 0],[0,-B/2,-D]),rad);


%SET MASS PROPERTIES
%assume body is box shape
setFrameMass(mdl,1,Mb,[0 0 Hb/2]',inertiaBox(Mb,L,Wb,Hb));

%assume wheels are cylinders
for i = mdl.wheelframeinds
    setFrameMass(mdl,i,Mw,[0 0 0]',inertiaCylinder(Mw,rad,Ww,2));
end



%FOR KINEMATIC MODEL
mdl.min_npic = nw;

mdl.dz_target = -.02*ones(1,nw);
mdl.tc_z = .1*ones(1,nw);
mdl.tc_j = .1;

%keep steer angles fixed when initializing terrain contact
steer_fi = namesToInds(mdl,{'s1','s2'});
for i = 1:length(steer_fi)
    set_isfixed(mdl,steer_fi(i),true);
end


% mdl.bsm_fh = @universal_bsm;
% mdl.bsm_p = zeros(6,1);

mdl.bsm_fh = [];

%FOR DYNAMIC MODEL

%set function handles

mdl.wgc_fh = @uniformWgc; %wheel-ground contact model
[~,~,fh]=feval(mdl.wgc_fh);

Kp=TotalMass*mdl.grav/(nw*-mdl.dz_target(1));

mdl.wgc_p = wgcParams(fh,Kp); %wheel force model parameters

%DEBUGGING
if strcmp(func2str(fh),'ishigamiWgc') || ...
   strcmp(func2str(fh),'ishigamiLUTWgc')
    
    disp([mfilename ': modifying wheel-ground contact parameters'])
    mdl.wgc_p(2) = rad;
    mdl.wgc_p(3) = Ww;
end

mdl.act_fh = @PIact;
mdl.act_p = [2e3 0 Inf]';



%ODE contact model parameters
[erp,cfm] = KpKdToErpCfm(Kp,Kp/20,.04);
fds = 1/(.1*Kp);

mdl.erp_z = erp*ones(1,nw);
mdl.cfm_z = cfm*ones(1,nw);
mdl.fds_x = fds*ones(1,nw);
mdl.fds_y = fds*ones(1,nw);
if ~opts.fix_front_axle_roll
    mdl.erp_j = .2;
    mdl.cfm_j = 1e-6;
end

%FOR BOTH
%set function handles
mdl.controller_fh = @zoeController;
% mdl.controller_fh = @zoeControllerSet;

if ~opts.fix_front_axle_roll
    mdl.hjc_fh = @zoeConstraints; %additional constraints, roll averaging
end

mdl.cov_p=zeros(4,1); %stochastic

%initialize state
orientation = [0 0 0]'*pi/180;
% orientation = eulerToQuat(orientation);

dlen=-1+length(orientation)+3; %ns-nf
ns = mdl.nf+dlen; %number of states
[isorient,ispos] = stateIs(ns);

state = zeros(ns,1);
state(isorient) = orientation;
state(ispos) = [0 0 rad+D-.01]';


% state(steer_fi+dlen) = 15*pi/180*[1 -1]; %DEBUGGING

%initialize velocity
% qvel = zeros(mdl.nf+5,1);

%DEBUGGING, nonzero initial velocity
[~,qvel]=feval(mdl.controller_fh,mdl,0,state);

% pathno = 1;
% [~,qvel]=feval(mdl.controller_fh,mdl,0,state,pathno);


%ANIMATION
cdir = fullfile(CADdir(),'Zoe');
% cdir = [CADdir() 'Zoe\']; %for windows

if nargout > 3
    anim = WmrAnimation();
    
%     fig_filename = fullfile('test','_autosave','anim_hgt_zoe');
    if 0 %load from .fig file (faster) 
%         loadHgt(anim,fig_filename)
    else
        addHgt(anim,[mdl.frames.parent_ind]);
        
        use_vrml = 1; %use VRML files
        if use_vrml 
            
            draw_edges = true;
            fix_lighting = true;
            alpha = 0;

            %Body
            i=1;
            makeVrmlHgt(anim.h_hgt(i),fullfile(cdir,'ZoeBody.wrl'),[],0,draw_edges,fix_lighting,alpha)

            %Front axle
            i = namesToInds(mdl,'s1');
            makeVrmlHgt(anim.h_hgt(i),fullfile(cdir,'ZoeFrontAxle.wrl'),[],0,draw_edges,fix_lighting,alpha)

            %Rear axle
            i = namesToInds(mdl,'s2'); 
            makeVrmlHgt(anim.h_hgt(i),fullfile(cdir,'ZoeRearAxle.wrl'),[],0,draw_edges,fix_lighting,alpha)

            %wheels
            i = namesToInds(mdl,'FL'); 
            makeVrmlHgt(anim.h_hgt(i),fullfile(cdir,'ZoeWheel.wrl'),[],0,draw_edges,fix_lighting,alpha)

            j = namesToInds(mdl,'BL');
            copyobj(anim.h_hgt(i),anim.h_hgt(j))
            
            j = namesToInds(mdl,'FR');
            copyobj(anim.h_hgt(i),anim.h_hgt(j))

            j = namesToInds(mdl,'BR');
            copyobj(anim.h_hgt(i),anim.h_hgt(j))
                     
        else

            %draw manually

            %Zoe body (octagon)

            bw = 1.15/2; %body width (half)
            bl = 1.70/2; %body length (half)
            dx = .3; %offsets from corners
            dy = .5;
            thk = .05; %axle thickness (half)

            pts(1,:) = [bl bl-dx -(bl-dx) -bl -bl -(bl-dx) bl-dx bl]; %X
            pts(2,:) = [-(bw-dy) -bw -bw -(bw-dy) bw-dy bw bw bw-dy]; %Y
            pts(:,end+1) = pts(:,1); %close the polygon

            C='none';
%             C=[0 1 1];
            makePolyExtrudeHgt(anim.h_hgt(1),pts,Hb,3,[0,0,Hb/2-thk]',1,1,C);

            %Zoe front axle
            i = namesToInds(mdl,'s1');

            dy = thk;
            aw = B/2-.01; %axle width (half)
            clear pts
            pts(1,:) = [-dy -(aw-thk) -aw -aw -(aw-thk) -dy dy (aw-thk) aw aw (aw-thk) dy]; %Y
            pts(2,:) = [thk -D+thk -D+thk -D-thk -D-thk -thk -thk -D-thk -D-thk -D+thk -D+thk thk]; %Z
            pts(:,end+1) = pts(:,1); %close the polygon

%             C=[0 1 1]/2;
            makePolyExtrudeHgt(anim.h_hgt(i),pts,thk*2,1,zeros(3,1),1,1,C);

            %Zoe rear axle
            j = namesToInds(mdl,'s2');

            copyobj(anim.h_hgt(i),anim.h_hgt(j))

            %Zoe wheels
            nv = 36; %num vertices
%             C = [1 1 1]/2;
            for i = mdl.wheelframeinds
                makeCircleHgt(anim.h_hgt(i),rad,2,zeros(3,1),nv,C);
            end
        end
    
        %total station prism
%         hprism=73*.0254;
%         C = 'none';
%         C = [1 1 1]/2;
%         makeCylinderHgt(anim.h_hgt(1),.02,hprism,3,[0 0 hprism/2]',36,true,C)

        %contact frames
        addHgt(anim,mdl.wheelframeinds)
        h = addQuiverZ(anim, mdl.nf+1:mdl.nf+nw, 0.75*rad);
        set(h,'LineWidth',1.5, 'Color','green')
    
        if use_vrml
%             saveHgt(anim,fig_filename)
        end
        
    end
    
    %DEBUGGING
    if 0
        %%
        h = addQuiverXYZ(anim, 1:mdl.nf+nw, .75*rad);
        set(h,'LineWidth',1.5,'Color','b')

        Names = {mdl.frames.name};
        cNames = Names(mdl.wheelframeinds);
        for i=1:length(cNames)
            cNames{i} = [cNames{i} '_c'];
        end
        Names = [Names cNames];

        h = addText(anim, 1:mdl.nf+nw, Names);
        set(h,'FontSize',16,'Color','m')
    end
    
    anim.ax_center = [0, 0, Hb/2]';
    d = 1.0*L;
    anim.ax_limits = [-d d -d d -d/2 d];
    
end


