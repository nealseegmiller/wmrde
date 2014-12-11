function [mdl, state, qvel, anim] = recbot()
%The RecBot is an instrumented John Deere Gator
%Specifications at:
%https://www.deere.com/wps/dcom/en_US/products/equipment/gator_utility_vehicles/traditional_utility_vehicles/t_series/te4x2electric/te4x2electric.page?

%min turn radius 7.3 m?

%conversion constants
intom = 2.54/100; %inch to meter

nw=4; %number of wheels
L = 1.770; %wheelbase (length)
B = [1.270 1.219]; %track width [front rear]
rad = [22 24]/2*intom; %wheel radius [front rear]

%Masses
TotalMass=667;
Mb = .95*TotalMass; %mass body
Mw = .05*TotalMass/nw; %mass wheel

%length, width, height of body
Lb=L;
Wb=min(B);
Hb=1.130-max(rad);

Ww = [9.5 12]*intom;

do_suspension = 0; %for dynamic sim only

%BUILD WmrModel OBJECT
mdl = WmrModel(nw+1+2+2*do_suspension);
addBodyFrame(mdl,'body');

%wheels:
zoff=rad(1)-rad(2)-.05; %difference between front/rear wheel radius, suspension deflection
if do_suspension
    addJointFrame(mdl,'suspL','body','PZ',0,poseToHT([0 0 0],[L, B(1)/2, 0]));
    addJointFrame(mdl,'suspR','body','PZ',0,poseToHT([0 0 0],[L, -B(1)/2, 0]));
    
    addJointFrame(mdl,'steerL','suspL','RZ',1,poseToHT([0 0 0],[0, 0, 0]));
    addJointFrame(mdl,'steerR','suspR','RZ',1,poseToHT([0 0 0],[0, 0, 0]));
else
    addJointFrame(mdl,'steerL','body','RZ',1,poseToHT([0 0 0],[L, B(1)/2, 0]));
    addJointFrame(mdl,'steerR','body','RZ',1,poseToHT([0 0 0],[L, -B(1)/2, 0]));
end

addWheelFrame(mdl,'FL','steerL',0,poseToHT([0 0 0],[0 0 zoff]),rad(1)); %front left
addWheelFrame(mdl,'FR','steerR',0,poseToHT([0 0 0],[0 0 zoff]),rad(1)); %front right

addWheelFrame(mdl,'BL','body',1,poseToHT([0 0 0],[0, B(2)/2, 0]),rad(2)); %rear left
addWheelFrame(mdl,'BR','body',1,poseToHT([0 0 0],[0, -B(2)/2, 0]),rad(2)); %rear right

%SET MASS PROPERTIES
%assume body is box shape
setFrameMass(mdl,1,Mb,[Lb/2 0 Hb/2]',inertiaBox(Mb,Lb,Wb,Hb));

%assume wheels are cylinders
wfi = mdl.wheelframeinds; %fl fr bl br
for i = wfi(1:2)
    setFrameMass(mdl,i,Mw,[0 0 0]',inertiaCylinder(Mw,rad(1),Ww(1),2));
end
for i = wfi(3:4)
    setFrameMass(mdl,i,Mw,[0 0 0]',inertiaCylinder(Mw,rad(2),Ww(2),2));
end

%FOR KINEMATIC MODEL
mdl.min_npic = nw;

mdl.dz_target = -.02*ones(1,nw);
mdl.tc_z = .1*ones(1,nw);
if do_suspension
    mdl.tc_j = .1*ones(1,2);
end

%keep front wheel rotations fixed when initializing terrain contact
fi = namesToInds(mdl,{'FL','FR'});
for i = 1:length(fi)
    set_isfixed(mdl,fi(i),true);
end

mdl.bsm_fh = @universalBsm;
mdl.bsm_p = zeros(6,1);

%FOR DYNAMIC MODEL
mdl.wgc_fh = @uniformWgc; %wheel-ground contact model

Kp = TotalMass*mdl.grav./(nw*-mdl.dz_target);

[~,~,fh]=feval(mdl.wgc_fh);
mdl.wgc_p = wgcParams(fh,Kp);

mdl.act_fh = @PI_act;
mdl.act_p = [250 0 Inf]';

%FOR BOTH
%set function handles
mdl.controller_fh = @ackermanController;

if do_suspension
    mdl.hjc_fh = @recbotConstraints;
end
mdl.cov_p = zeros(4,1);


%initialize state
orientation = [0 0 0]'*pi/180;
% orientation = euler_to_quat(orientation);

dlen=-1+length(orientation)+3; %ns-nf
ns = mdl.nf+dlen; %number of states
[isorient,ispos] = stateIs(ns);

state = zeros(ns,1);
state(isorient) = orientation;
state(ispos) = [-5 0 max(rad)]';

%initialize velocity
% qvel = zeros(mdl.nf+5,1);

[~,qvel]=feval(mdl.controller_fh,mdl,0,state); %nonzero initial velocity



cdir = [CADdir 'RecBot\'];

if nargout > 3
    anim = WmrAnimation();
    
    if 0 %load from .fig file (faster) 
        loadHgt(anim,'_autosave/anim_hgt_gator')
    else
        addHgt(anim,mdl.parent_ind);
        
        use_vrml = 0; %use VRML files
        if use_vrml 
            
            %construct from VRML files
            showedges = 0;
            fixlighting = 1;
            alpha = 1;
            alphaw = 1;

            %Body
            i=1;
            makeVrmlHgt(anim.h_hgt(i),[cdir 'RecBotBody.wrl'],[],0,showedges,fixlighting,alpha)

            %wheels
            i = namesToInds(mdl,'FL');
            makeVrmlHgt(anim.h_hgt(i),[cdir 'RecBotWheelFront.wrl'],[],0,showedges,fixlighting,alphaw)

            j = namesToInds(mdl,'FR');
            copyobj(anim.h_hgt(i),anim.h_hgt(j))

            i = namesToInds(mdl,'BL');
            makeVrmlHgt(anim.h_hgt(i),[cdir 'RecBotWheelRear.wrl'],[],0,showedges,fixlighting,alphaw)

            j = namesToInds(mdl,'BR');
            copyobj(anim.h_hgt(i),anim.h_hgt(j))
        
        else %draw manually

            %vehicle body, a box
            C = 'none';
%             C = [1 1 1]/5;
            makeBoxHgt(anim.h_hgt(1),Lb,Wb,Hb,[Lb/2,0,Hb/2]',C);

            %wheels, circles
            nv = 36; %num vertices
%             C = [1 1 1]/2;
            for i = wfi(1:2)
                makeCircleHgt(anim.h_hgt(i),rad(1),2,zeros(3,1),nv,C);
            end
            for i = wfi(3:4)
                makeCircleHgt(anim.h_hgt(i),rad(2),2,zeros(3,1),nv,C);
            end
        end
        
        %contact frames
        addHgt(anim,mdl.wheelframeinds)
        h = addQuiverZ(anim, mdl.nf+1:mdl.nf+nw, 0.75*min(rad));
        set(h,'LineWidth',1.5, 'Color','green')
    
        if use_vrml
            saveHgt(anim,'_autosave/anim_hgt_gator')
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
    
    anim.ax_center = [Lb/2 0 Hb/2]';
    d = L;
    anim.ax_limits = [-d d -d d -d/2 d/2];
end

% blah=1;






