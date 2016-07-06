function [mdl, state, qvel, anim] = talon() 

%dimensions from:
%https://www.qinetiq-na.com/wp-content/uploads/brochure_talon.pdf


%conversion constants
intom = 2.54/100; %inch to meter
lbtokg = .453592; %pounds to kg

%???
rad = 5*intom; %sprocket radius
Wt = 6*intom; %width of track

%DEBUGGING
rad1 = 1.0*rad;
rad2 = rad;

nt=2; %number of tracks
L = (34*intom)-rad1-rad2; %distance between forward/rear sprockets
B = (22.5*intom)-Wt; %distance between left/right track centers


%Masses
TotalMass = 115*lbtokg;
Mb = TotalMass; %mass of body
Ms = 1e-6; %mass of sprocket, epsilon
Iy = .05*TotalMass*rad^2/2; %rotational inertia of track about y axis

%length, width, height of body
Lb = L;
Wb = B - Wt;
Hb = 1.5*rad;


%BUILD KINEMATIC TREE
mdl = WmrModel(nt+1);
addBodyFrame(mdl,'body');

%for tracks


addSprocketFrame(mdl,'TL','body',1,poseToHT([0 0 0],[-L/2,  B/2, 0]),rad1,rad2,L); %left
addSprocketFrame(mdl,'TR','body',1,poseToHT([0 0 0],[-L/2, -B/2, 0]),rad1,rad2,L); %right


%SET MASS PROPERTIES
%assume body is box shape
setFrameMass(mdl,1,Mb,[0 0 0]',inertiaBox(Mb,Lb,Wb,Hb));

%sprocket
for i = mdl.sprocketframeinds
    setFrameMass(mdl,i,Ms,[0 0 0]',diag([0 Iy 0]) );
end

%FOR KINEMATIC MODEL
mdl.update_Is();
mdl.min_npic = 3;

mdl.dz_target = -.01*ones(1,nt);
mdl.tc_z = .1*ones(1,nt);

mdl.bsm_fh = [];
mdl.bsm_p = [];

% mdl.bsm_fh = @universalBsm;
% mdl.bsm_p = zeros(6,1);


%FOR DYNAMIC MODEL
mdl.wgc_fh = @uniformWgc; %wheel slip model
[~,~,fh]=feval(mdl.wgc_fh);

npflat = 3; %number of points on flat section of track
Kp = TotalMass*mdl.grav/(nt*npflat*-mdl.dz_target(1));

mdl.wgc_p = wgcParams(fh,Kp); %wheel-ground contact model parameters

mdl.act_fh = @PIact;
mdl.act_p = [2e3 0 Inf]';

%ODE contact model parameters
[erp,cfm] = KpKdToErpCfm(Kp,Kp/20,.04);
fds = 1/(.1*Kp);

mdl.erp_z = erp*ones(1,nt);
mdl.cfm_z = cfm*ones(1,nt);
mdl.fds_x = fds*ones(1,nt);
mdl.fds_y = fds*ones(1,nt);
mdl.erp_j = .2;
mdl.cfm_j = 1e-6;

%FOR BOTH
mdl.controller_fh = @skidsteerController;
mdl.cov_p = zeros(4,1);


%initialize state
orientation = [0 0 0]'*pi/180;
position = [0 -.9 (rad-.01)]';

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

cdir = fullfile(CADdir(),'Talon');

if nargout > 3
    anim = WmrAnimation();
    
%     fig_filename = fullfile('test','_autosave','anim_hgt_talon');
    if 0 %load from .fig file (faster) 
%         loadHgt(anim,fig_filename)
    else
        addHgt(anim,[mdl.frames.parent_ind]);
        
        use_vrml = 1; %use VRML files
        if use_vrml

            %construct from VRML files
            draw_edges = 1;
            fix_lighting = 1;
            
            alpha = 0;

            %Body
            i=1;
            makeVrmlHgt(anim.h_hgt(i),fullfile(cdir,'TalonBody.wrl'),[],0,draw_edges,fix_lighting,alpha)

            %wheels
            for i = mdl.sprocketframeinds
                makeVrmlHgt(anim.h_hgt(i),fullfile(cdir,'TalonSprocket.wrl'),[],0,draw_edges,fix_lighting,alpha)
            end

        else %draw manually

            %construct from scratch
            %vehicle body, a box
            C = 'none';
%             C = [1 1 1]/5;
            makeBoxHgt(anim.h_hgt(1),Lb,Wb,Hb,[0 0 0]',C);

            %sprockets, circles
            nv = 36; %num vertices

            for i = mdl.sprocketframeinds;
                makeCircleHgt(anim.h_hgt(i),rad1,2,zeros(3,1),nv,C);
            end

            
        end
        
        %track frames
        contacts = initTrackContactGeom(mdl.frames(mdl.sprocketframeinds));
        for tno = 1:nt
            fi = mdl.sprocketframeinds(tno);
            parent_fi = mdl.frames(fi).parent_ind;
            addHgt(anim,parent_fi);
            i = mdl.nf + tno; %hgtransform handle index
            
            updateHgt(anim,mdl.frames(fi).HT_parent_jd0,i); %doesn't change

            %tracks
            if use_vrml
                makeVrmlHgt(anim.h_hgt(i),fullfile(cdir,'TalonTrack.wrl'),[],0,draw_edges,fix_lighting,alpha)
            else
                makeTrackHgt(anim.h_hgt(i),rad1,rad2,L,Wt,zeros(3,1),true,C)
            end
            
            %contact points
            np = size(contacts(tno).HT_track,3);
            cp = zeros(3,np);
            cp(:) = contacts(tno).HT_track(1:3,4,:);
            h = plot3(cp(1,:), cp(2,:), cp(3,:), 'k.');
            set(h,'Parent',anim.h_hgt(i));
            
            %normal vectors
            anim.h_quiver(tno) = drawQuivers(zeros(3,1),zeros(3,1));
            set(anim.h_quiver(tno),'Parent',anim.h_hgt(i),'Color','g');
        end
        

        if use_vrml
%             saveHgt(anim,fig_filename)
        end


    end
    
    anim.ax_center = [0 0 0]';
    d = 1.5*L;
    anim.ax_limits = [-d d -d d -d/2 d/2];
    

end

