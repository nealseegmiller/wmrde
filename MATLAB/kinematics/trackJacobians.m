function A = trackJacobians(mdl,HT_world,contacts)
%compute track Jacobians. A*qvel = v
%where v comprises contact point velocities (relative to world, in contact coords)
%INPUTS
%mdl:           WmrModel object
%HT_world:      4x4xnf matrix, transforms to world coords for WmrModel frames
%contacts:      1 x nt array of TrackContactGeom objects

%get from WmrModel
nt = mdl.nt;

nf = mdl.nf;
nv = nf-1+6; %length of qvel
parent_ind = [mdl.frames.parent_ind];
dof_type = [mdl.frames.dof_type];

R_body_to_world = HT_world(1:3,1:3,1);

Atrack = cell(1,nt); %Jacobians for each track

%for indexing spatial vectors
ang = 1:3;
lin = 4:6;

if 1
    %C++ version, similar to wheelJacobians
    
    for tno = 1:nt %track index
        
        %for this track:
        incontactinds = find(contacts(tno).incontact);
        npic = length(incontactinds); %number of points in contact
        
        Atrack{tno} = zeros(3*npic,nv); 
        
        if npic == 0
            continue
        end

        sprocket_fi = mdl.sprocketframeinds(tno); %sprocket frame index
        HT_track_to_world = HT_world(:,:, parent_ind(sprocket_fi)) * mdl.frames(sprocket_fi).HT_parent_jd0;
        rad = mdl.frames(sprocket_fi).rad;

        for i = 1:npic
            pno = incontactinds(i);
            
            Apt = zeros(3,nv); %Jacobian for a single point
            
            HT_contact_to_world = HT_track_to_world * contacts(tno).HT_track(:,:,pno);
            
            %sprocket angular rate
            Apt(:,sprocket_fi+5) = -HT_contact_to_world(1:3,1)*rad;
            
            %traverse up the chain
            fi = parent_ind(sprocket_fi);
            
            while fi > 1 %stop at body frame, handled below
                vi = fi+5; %index in qvel
                if any(dof_type(fi) == [1 2 3])
                    %revolute joint
                    %cross rotation axis with translation vector
                    ax = HT_world(1:3,dof_type(fi),fi);
                    r = HT_contact_to_world(1:3,4) - HT_world(1:3,4,fi);
                    Apt(:,vi) = cross3(ax,r); 
                else
                    %prismatic joint
                    ax = HT_world(1:3,dof_type(fi)-3,fi);
                    Apt(:,vi) = ax;
                end
                fi = parent_ind(fi); %update to next parent
            end

            %body joint
            %right multiply by R_body_to_world because body frame velocity is in body coords in qvel
            Apt(:,lin) = R_body_to_world;
            r = HT_contact_to_world(1:3,4) - HT_world(1:3,4,1);
            Apt(:,ang) = skew3(r)'*R_body_to_world;
            
            %rotate into track frame coords
            Apt = HT_contact_to_world(1:3,1:3)'*Apt;
            
            I = (i-1)*3 + (1:3); %row indices
            Atrack{tno}(I,:) = Apt;

        end
        
    end
    
    
else
    
    %MATLAB version, fewer for loops
    %requires precomputed .dvc_dvt
    
    for tno = 1:nt %track index

        incontact = contacts(tno).incontact;

        npic = sum(incontact);
        Atrack{tno} = zeros(3*npic,nv);
        
        if npic == 0
            continue
        end

        %Jacobian of [track frame spatial vel; sprocket angular rate] wrt qvel
        dvt_dqvel = zeros(7,nv);

        fi = mdl.sprocketframeinds(tno); %sprocket frame index
        HT_track_to_world = HT_world(:,:,parent_ind(fi)) * mdl.frames(fi).HT_parent_jd0;
        dvt_dqvel(7,fi+5) = 1; %sprocket angular rate

        fi = parent_ind(fi);

        %traverse up the chain
        while fi > 1 %stop at body frame, handled below
            vi = fi+5; %index in qvel
            if any(dof_type(fi) == [1 2 3])
                %revolute joint
                %cross rotation axis with translation vector
                ax = HT_world(1:3,dof_type(fi),fi);
                dvt_dqvel(ang,vi) = ax;
                r = HT_track_to_world(1:3,4) - HT_world(1:3,4,fi);
                dvt_dqvel(lin,vi) = cross3(ax,r); 
            else
                %prismatic joint
                ax = HT_world(1:3,dof_type(fi)-3,fi);
                dvt_dqvel(lin,vi) = ax;
            end
            fi = parent_ind(fi); %update to next parent
        end

        %body joint
        %right multiply by R_body_to_world because body frame velocity is in body coords in qvel
        dvt_dqvel(ang,ang) = R_body_to_world;
        r = HT_track_to_world(1:3,4) - HT_world(1:3,4,1);
        dvt_dqvel(lin,ang) = skew3(r)'*R_body_to_world;
        dvt_dqvel(lin,lin) = R_body_to_world;

        %rotate into track frame coords
        R_world_to_track = HT_track_to_world(1:3,1:3)';
        dvt_dqvel(ang,:) = R_world_to_track * dvt_dqvel(ang,:);
        dvt_dqvel(lin,:) = R_world_to_track * dvt_dqvel(lin,:);

        %row indices for points in contact
        np = size(contacts(tno).HT_track,3);
        I = reshape(1:3*np,3,np);
        I = I(:,incontact);

        Atrack{tno} = contacts(tno).dvc_dvt(I,:) * dvt_dqvel;
    end

end

A = cat(1,Atrack{:});

