function [mdl, state, qvel, anim] = landtamer()

%obtained dimensions from diagram at:
%http://www.landtamer.com/product/1/95
%http://www.remoteaccessvehicles.com/graphics/LT-brochure-ExtremeLogistics.pdf

%conversion constants
intom = 2.54/100; %inch to meter
lbtokg = .453592; %pounds to kg

nw=6; %number of wheels
L = 84*intom; %wheelbase (length)
B = (72.5-8)*intom; %track (width)
rad = 33/2*intom; %wheel radius

%Masses
TotalMass=3225*lbtokg;
Mb = .95*TotalMass; %mass body
Mw = .05*TotalMass/nw; %mass wheel

%length, width, height of body
Lb=L;
Wb=B;
Hb=32.6*intom;

Ww=10*intom; %width of wheel

%BUILD KINEMATIC TREE
mdl = WmrModel(nw+1);
addBodyFrame(mdl,'body');

%wheels:
addWheelFrame(mdl,'FL','body',1,poseToHT([0 0 0],[L/2,  B/2, 0]),rad); %front left
addWheelFrame(mdl,'FR','body',1,poseToHT([0 0 0],[L/2, -B/2, 0]),rad); %front right

addWheelFrame(mdl,'ML','body',1,poseToHT([0 0 0],[0,    B/2, 0]),rad); %mid left
addWheelFrame(mdl,'MR','body',1,poseToHT([0 0 0],[0,   -B/2, 0]),rad); %mid right

addWheelFrame(mdl,'BL','body',1,poseToHT([0 0 0],[-L/2, B/2, 0]),rad); %rear left
addWheelFrame(mdl,'BR','body',1,poseToHT([0 0 0],[-L/2,-B/2, 0]),rad); %rear right

%SET MASS PROPERTIES
%assume body is box shape
setFrameMass(mdl,1,Mb,[0 0 Hb/2]',inertiaBox(Mb,Lb,Wb,Hb));
% setFrameMass(mdl,1,Mb,[0 0 Hb]',inertiaBox(Mb,Lb,Wb,Hb)); disp([mfilename ': DEBUGGING, cg offset'])

%assume wheels are cylinders
for i = mdl.wheelframeinds
    setFrameMass(mdl,i,Mw,[0 0 0]',inertiaCylinder(Mw,rad,Ww,2));
end


%FOR KINEMATIC MODEL
mdl.min_npic = 3;
mdl.update_Is();

mdl.dz_target = -.02*ones(1,nw);
% mdl.dz_target = -.05*ones(1,nw); %?
mdl.tc_z = .1*ones(1,nw);

mdl.bsm_fh = @universalBsm;
mdl.bsm_p = zeros(6,1);

%FOR DYNAMIC MODEL

mdl.wgc_fh = @uniformWgc; %wheel slip model
[~,~,fh]=feval(mdl.wgc_fh);

Kp = TotalMass*mdl.grav/(nw*-mdl.dz_target(1));

mdl.wgc_p = wgcParams(fh,Kp); %wheel-ground contact model parameters

mdl.act_fh = @PIact;
mdl.act_p = [2e3 0 Inf]';

%FOR BOTH
mdl.controller_fh = @skidsteerController;
mdl.cov_p = zeros(4,1);

%initialize state
orientation = [0 0 0]'*pi/180; %Euler angles
% position = [0 0 rad]';
position = [-12 -7 rad]';


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

[~,qvel]=feval(mdl.controller_fh,mdl,0,state); %nonzero initial velocity

%ANIMATION
cdir = [CADdir() 'LandTamer/'];

if nargout > 3
    anim = WmrAnimation();
    
    if 1 %load from .fig file (faster) 
        loadHgt(anim,'_autosave/anim_hgt_landtamer')
    else
        addHgt(anim,[mdl.frames.parent_ind]);
        
        use_vrml = 1; %use VRML files
        if use_vrml

            %construct from VRML files
            draw_edges = 0;
            fix_lighting = 1;
            
            alpha = 1;
            alpha_wheel = 1;

            %Body
            i=1;
            makeVrmlHgt(anim.h_hgt(i),[cdir 'LandTamerBody.wrl'],[],0,draw_edges,fix_lighting,alpha)

            %wheels
            i = namesToInds(mdl,'FL');
            makeVrmlHgt(anim.h_hgt(i),[cdir 'LandTamerWheel.wrl'],[],0,draw_edges,fix_lighting,alpha_wheel)

            j = namesToInds(mdl,'ML');
            copyobj(anim.h_hgt(i),anim.h_hgt(j))

            j = namesToInds(mdl,'BL');
            copyobj(anim.h_hgt(i),anim.h_hgt(j))

            j = namesToInds(mdl,'FR');
            copyobj(anim.h_hgt(i),anim.h_hgt(j))

            j = namesToInds(mdl,'MR');
            copyobj(anim.h_hgt(i),anim.h_hgt(j))

            j = namesToInds(mdl,'BR');
            copyobj(anim.h_hgt(i),anim.h_hgt(j))


        else %draw manually
        
        

            %construct from scratch
            %vehicle body, a box
            C = 'none';
%             C = [1 1 1]/5;
            makeBoxHgt(anim.h_hgt(1),Lb,Wb,Hb,[0,0,Hb/2]',C);

            %wheels, circles
            nv = 36; %num vertices
%             C = [1 1 1]/2;
            for i = mdl.wheelframeinds
                makeCircleHgt(anim.h_hgt(i),rad,2,zeros(3,1),nv,C);
            end
            
        end
        
        %contact frames
        addHgt(anim,mdl.wheelframeinds)
        h = addQuiverZ(anim, mdl.nf+1:mdl.nf+nw, 0.75*rad);
        set(h,'LineWidth',1.5, 'Color','green')

        if use_vrml
            saveHgt(anim,'_autosave/anim_hgt_landtamer')
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
    
    anim.ax_center = [0 0 Hb/2]';
    d = 1.5*L;
    anim.ax_limits = [-d d -d d -d/2 d/2];
    

end





