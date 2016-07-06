function [mdl, state, qvel, anim] = rocky()

%sources
%Balaram, Ohm, and Ivlev, The Rocky 7 Mars rover prototype, IROS 1996
%Tarokh and McDermott, Kinematics modeling and analyses of articulated rovers, IEEE Transactions on Robotics 2005

%Dimensions
k1 = 10.5/100;      %vertical offset between (R)over reference to (D)ifferential
k2 = 12.075/100;    %forward offset between R and D
k3 = 20/100;        %horizontal distance between D and wheels
k4 = 28.8/100;      %distance from D to steering axis of front wheels
k5 = 12.5/100;      %height of D from wheel axles
k6 = 16/100;        %length of link from rocker joint to bogie joint
k7 = 6.9/100;       %length from bogie joint to front/rear bogie
k8 = 2/100;         %height of bogie joint from wheel axles
k9 = 139;           %angle of link between rocker and bogie joints (degrees)
k10 = 6.5/100;      %wheel radius

k6x = sind(k9-90)*k6;
k6z = cosd(k9-90)*k6;



%Mass properties
Mb = 11; %mass of body
Mw = .5; %mass of each wheel

%length, width, height of body
Lb = .40;
Wb = .26;
Hb = .16;
rad = k10; %wheel radius
Ww = .08; %wheel width

%body frame center of mass (in body coordinates)
%(from CAD)
cmx = 18/100;
cmz = 15.5/100;


%BUILD KINEMATIC TREE
mdl = WmrModel(13); %robot reference/body frame
addBodyFrame(mdl,'R');

%left side
addJointFrame(mdl,'D1','R','RY',0,poseToHT([0 0 0],[k2, k3, k1])); %differential frame, left

addJointFrame(mdl,'S1','D1','RZ',1,poseToHT([0 0 0],[k4, 0, 0])); %steering frame, front left
addWheelFrame(mdl,'A1','S1',1,poseToHT([0 0 0],[0, 0, -k5]),rad); %wheel frame, front left

addJointFrame(mdl,'B1','D1','RY',0,poseToHT([0 0 0],[-k6x, 0, -k6z])); %bogie frame, left
addWheelFrame(mdl,'A3','B1',1,poseToHT([0 0 0],[ k7, 0, -k8]),rad); %wheel frame, left bogie-front
addWheelFrame(mdl,'A5','B1',1,poseToHT([0 0 0],[-k7, 0, -k8]),rad); %wheel frame, left bogie-rear

%right side
addJointFrame(mdl,'D2','R','RY',0,poseToHT([0 0 0],[k2, -k3, k1])); %differential frame, right

addJointFrame(mdl,'S2','D2','RZ',1,poseToHT([0 0 0],[k4, 0, 0])); %steering frame, front right
addWheelFrame(mdl,'A2','S2',1,poseToHT([0 0 0],[0, 0, -k5]),rad); %wheel frame, front right

addJointFrame(mdl,'B2','D2','RY',0,poseToHT([0 0 0],[-k6x, 0, -k6z])); %bogie frame, right
addWheelFrame(mdl,'A4','B2',1,poseToHT([0 0 0],[ k7, 0, -k8]),rad); %wheel frame, right bogie-front
addWheelFrame(mdl,'A6','B2',1,poseToHT([0 0 0],[-k7, 0, -k8]),rad); %wheel frame, right bogie-rear


%SET MASS PROPERTIES
nw = mdl.nw;
TotalMass = Mb + nw*Mw;

%assume body is box shape
setFrameMass(mdl,1,Mb,[cmx 0 cmz]',inertiaBox(Mb,Lb,Wb,Hb));

%assume wheels are cylinders
for i = mdl.wheelframeinds
    setFrameMass(mdl,i,Mw,[0 0 0]',inertiaCylinder(Mw,rad,Ww,2));
end


%FOR KINEMATIC MODEL

%all wheels must be in contact
mdl.min_npic = nw;

% mdl.dz_target = -.002*ones(1,nw);
mdl.dz_target = -.02*ones(1,nw);
mdl.tc_z = .1*ones(1,nw);
mdl.tc_j = .1;


%FOR DYNAMIC MODEL

mdl.wgc_fh = @uniformWgc; %wheel-ground contact model
[~,~,fh]=feval(mdl.wgc_fh);
Kp = TotalMass*mdl.grav/(nw*-mdl.dz_target(1));
mdl.wgc_p = wgcParams(fh,Kp); %wheel force model parameters

%DEBUGGING
if strcmp(func2str(fh),'ishigamiWgc') || ...
   strcmp(func2str(fh),'ishigamiLUTWgc')
    
    disp([mfilename ': modifying wheel-ground contact parameters'])
    mdl.wgc_p(2) = rad;
    mdl.wgc_p(3) = Ww;
end

mdl.act_fh = @PIact;
mdl.act_p = [5 0 Inf]';

%ODE contact model parameters
[erp,cfm] = KpKdToErpCfm(Kp,Kp/20,.04);
fds = 1/(.1*Kp);

mdl.erp_z = erp*ones(1,nw);
mdl.cfm_z = cfm*ones(1,nw);
mdl.fds_x = fds*ones(1,nw);
mdl.fds_y = fds*ones(1,nw);
mdl.erp_j = .2;
mdl.cfm_j = 1e-6;

%FOR BOTH
%set function handles
mdl.controller_fh = @rockyController;
mdl.hjc_fh = @rockyConstraints;

%initialize state
orientation = [0 0 0]'*pi/180; %Euler angles
position = [1.1 0 (k10+k8+k6z-k1)]';

olen = SIZEORIENT();
dlen = -1+olen+3;
ns = mdl.nf + dlen; %number of states
[isorient,ispos] = stateIs(ns);

if olen == 4
    %TODO
end

state = zeros(ns,1);
state(isorient) = orientation;
state(ispos) = position;

%initialize velocity
% qvel = zeros(mdl.nf+5,1);
[~,qvel]=feval(mdl.controller_fh,mdl,0,state); %DEBUGGING, nonzero initial velocity



%ANIMATION
cdir = fullfile(CADdir(),'Rocky7');

if nargout > 3
    anim = WmrAnimation();
    
%     fig_filename = fullfile('test','_autosave','anim_hgt_rocky');
    if 0 %load from .fig file (faster) 
%         loadHgt(anim,fig_filename)
    else
        addHgt(anim,[mdl.frames.parent_ind]);
        
        use_vrml = 1; %use VRML files
        if use_vrml
            
            draw_edges = 1;
            fix_lighting = 1;
            alpha = 0;
            alpha_body = 0;
            alpha_wheel = 0;
            HT = [];

            %Body
            makeVrmlHgt(anim.h_hgt(1),fullfile(cdir,'Rocky7Body.wrl'),HT,0,draw_edges,fix_lighting,alpha_body);

            %Left Rocker
            i = namesToInds(mdl,'D1');
            makeVrmlHgt(anim.h_hgt(i),fullfile(cdir,'Rocky7Rocker.wrl'),HT,0,draw_edges,fix_lighting,alpha);

            %Right Rocker
            i = namesToInds(mdl,'D2');
            makeVrmlHgt(anim.h_hgt(i),fullfile(cdir,'Rocky7Rocker.wrl'),HT,2,draw_edges,fix_lighting,alpha);

            %Left Bogie
            i = namesToInds(mdl,'B1');
            makeVrmlHgt(anim.h_hgt(i),fullfile(cdir,'Rocky7Bogie.wrl'),HT,0,draw_edges,fix_lighting,alpha);

            %Right Bogie
            i = namesToInds(mdl,'B2');
            makeVrmlHgt(anim.h_hgt(i),fullfile(cdir,'Rocky7Bogie.wrl'),HT,2,draw_edges,fix_lighting,alpha);

            %Left Steering Bracket
            i = namesToInds(mdl,'S1');
            makeVrmlHgt(anim.h_hgt(i),fullfile(cdir,'Rocky7Bracket.wrl'),HT,0,draw_edges,fix_lighting,alpha);

            %Left Wheels
            i = namesToInds(mdl,'A1');
            makeVrmlHgt(anim.h_hgt(i),fullfile(cdir,'Rocky7Wheel.wrl'),HT,0,draw_edges,fix_lighting,alpha_wheel);

            for j = namesToInds(mdl,{'A3','A5'});
                copyobj(anim.h_hgt(i),anim.h_hgt(j))
            end

            %Right Steering Bracket
            i = namesToInds(mdl,'S2');
            makeVrmlHgt(anim.h_hgt(i),fullfile(cdir,'Rocky7Bracket.wrl'),HT,2,draw_edges,fix_lighting,alpha);

            %Right Wheels
            i = namesToInds(mdl,'A2');
            makeVrmlHgt(anim.h_hgt(i),fullfile(cdir,'Rocky7Wheel.wrl'),HT,2,draw_edges,fix_lighting,alpha_wheel);

            for j = namesToInds(mdl,{'A4','A6'});
                copyobj(anim.h_hgt(i),anim.h_hgt(j))
            end      
            
        else %draw manually
        
            %Rocky body, a box
            C = 'none';
%             C = [0 1 1];
            makeBoxHgt(anim.h_hgt(1),Lb,Wb,Hb,[cmx 0 cmz]',C);

            yoff=5.4/100;
            thk=2/100;
            %Rocker
            i = namesToInds(mdl,'D1');
            pts(1,:) = [k4, 0, -k6x];
            pts(2,:) = [0, 0, -k6z];
%             C = [0 1 1]/2;
            makeLineExtrudeHgt(anim.h_hgt(i),pts,thk,thk,2,[0,-yoff,0]',1,C)

            %2nd Rocker
            i = namesToInds(mdl,'D2');
            makeLineExtrudeHgt(anim.h_hgt(i),pts,thk,thk,2,[0,yoff,0]',1,C)

            %Bogie
            i = namesToInds(mdl,'B1');
            pts(1,:) = [k7, 0, -k7];
            pts(2,:) = [-k8, 0, -k8];
            makeLineExtrudeHgt(anim.h_hgt(i),pts,thk,thk,2,[0,-yoff,0]',1,C)

            %2nd Bogie
            i = namesToInds(mdl,'B2');
            makeLineExtrudeHgt(anim.h_hgt(i),pts,thk,thk,2,[0,yoff,0]',1,C)

            %Wheels
            nv = 36; %num vertices
%             C = [1 1 1]/2;
            for i = mdl.wheelframeinds
                makeCircleHgt(anim.h_hgt(i),rad,2,zeros(3,1),nv,C);
%                 makeCylinderHgt(anim.h_hgt(i),rad,Ww,2,zeros(3,1),nv,true,C)
            end
            
        end
        
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
    
    anim.ax_center = [cmx 0 cmz]';
    d = .5;
    anim.ax_limits = [-d d -d d -d d];
end







