function H = jointSpaceInertia(mdl,Xup,Is_subt)
%Compute joint space bias force per Featherstone Composite Rigid Body Algorithm
%adapted from HandC.m in spatial_v2 library
%http://users.cecs.anu.edu.au/~roy/spatial/
%INPUT
%obj:       WmrModel object
%Xup:       6x6x nf, Plucker transforms, Xup(:,:,i) transforms spatial motion vector from frame parent(i) to i coords
%Is_subt:   6x6x nf, spatial inertia of subtree rooted at each frame
%OUTPUT
%H:         nv x nv, joint space inertia matrix

%get from WmrModel object
nf = mdl.nf;
nv = mdl.nf+5;
dof_type = [mdl.frames.dof_type];
parent_ind = [mdl.frames.parent_ind];

H = zeros(nv); %Generalized inertia matrix

%COMPUTE H
H(1:6,1:6) = Is_subt(:,:,1); %body frame
for fi = 2:nf %frame index
    vi = fi+5; %vel index (in qvel)

    force_h = Is_subt(:,dof_type(fi),fi);
    H(vi,vi) = force_h(dof_type(fi));

    fj = fi; %2nd frame index
    while fj > 1
        force_h = Xup(:,:,fj)' * force_h;
        fj = parent_ind(fj);

        if fj == 1
            %body frame 'joint'
            vj = 1:6; %2nd vel index
            H(vj,vi) = force_h;
        else
            vj = fj+5; %2nd vel index
            H(vj,vi) = force_h(dof_type(fj));
        end
        H(vi,vj) = H(vj,vi)'; %symmetry
    end
end
H = (H+H')/2; %make symmetric


