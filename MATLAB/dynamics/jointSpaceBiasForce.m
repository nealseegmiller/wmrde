function C = jointSpaceBiasForce(mdl, Xup, qvel)
%Compute joint space bias force per Featherstone Recursive Newton Euler Algorithm
%adapted from HandC.m in spatial_v2 library
%http://users.cecs.anu.edu.au/~roy/spatial/
%INPUT
%mdl:       WmrModel object
%Xup:       6x6x nf, Plucker transforms, Xup(:,:,i) transforms spatial motion vector from frame parent(i) to i coords
%qvel:      nv x 1, joint space velocity
%OUTPUT
%C:         nv x 1, joint space bias force

%get from WmrModel object
nf = mdl.nf;
nv = nf-1+6;

dof_type = [mdl.frames.dof_type];
parent_ind = [mdl.frames.parent_ind];

% mass = [mdl.frames.mass];
Is = cat(3,mdl.frames.Is);

%init outputs
C = zeros(nv,1); %joint-space bias force

%variable names from HandC.m in spatial_v2 library
%spatial vectors for each frame, in link coords:
v = zeros(6,nf); %velocity
avp = zeros(6,nf); %acceleration (includes Coriolis, centripetal, and gravity terms)
fvp = zeros(6,nf); %force

%spatial acc of gravity in body coords
a_grav = zeros(6,1);
a_grav(4:6) = Xup(1:3,3,1)*-mdl.grav; %world z axis in body coords * grav acc


%body frame
v(:,1) = qvel(1:6); %spatial velocity of body frame
avp(:,1) = -a_grav; %accelerating the base is faster than explicitly adding gravity force to each body

for fi = 1:nf

    if fi > 1
        pfi = parent_ind(fi); %parent frame ind
        vJ = zeros(6,1); %spatial joint velocity
        vJ(dof_type(fi)) = qvel(fi+5); 
        
        v(:,fi) = Xup(:,:,fi)*v(:,pfi) + vJ;
        avp(:,fi) = Xup(:,:,fi)*avp(:,pfi) + crossMotion(v(:,fi))*vJ;  
    end

%     if mass(fi) > 0
        fvp(:,fi) = Is(:,:,fi)*avp(:,fi) + crossForce(v(:,fi))*Is(:,:,fi)*v(:,fi);
%     end
end


for fi = nf:-1:2
    C(fi+5) = fvp(dof_type(fi),fi);
    
    pfi = parent_ind(fi);
    %Transform force to parent frame coords
    fvp(:,pfi) = fvp(:,pfi) + Xup(:,:,fi)'*fvp(:,fi); 
end
C(1:6) = fvp(:,1); %body frame

% return





