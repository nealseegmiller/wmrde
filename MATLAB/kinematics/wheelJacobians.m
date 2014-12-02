function A = wheelJacobians(mdl,HT_world,contacts)
%compute wheel Jacobians. A*qvel = v
%where v comprises contact point velocities (relative to world, in contact coords)
%INPUTS
%mdl:           WmrModel object
%HT_world:      4x4xnf matrix, transforms to world coords for WmrModel frames
%contacts:      1 x nw array of WheelContactGeom objects

%get from WmrModel
nf = mdl.nf;
nv = nf-1+6; %length of qvel
parent_ind = [mdl.frames.parent_ind];
dof_type = [mdl.frames.dof_type];

R_body_to_world = HT_world(1:3,1:3,1);

incontact = [contacts.incontact];
A = zeros(sum(incontact)*3,nv);

%for indexing spatial vectors
ang = 1:3;
lin = 4:6;

i = 0; %counter for wheels in contact
for wno = find(incontact) %wheel index
    
    HT_contact_to_world = contacts(wno).HT_world;
    
    Apt = zeros(3,nv); %Jacobian for a single point
        
    fi = mdl.wheelframeinds(wno); %index of parent wheel
    %traverse up the chain
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
    r = HT_contact_to_world(1:3,4) - HT_world(1:3,4,1);
    Apt(:,ang) = skew3(r)'*R_body_to_world;
    Apt(:,lin) = R_body_to_world;
    
    %rotate into contact frame coords
    Apt = HT_contact_to_world(1:3,1:3)'*Apt;

    I = i*3 + (1:3); %row indices
    A(I,:) = Apt;
    
    i = i + 1;
end


