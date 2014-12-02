function [HT_parent, HT_world] = stateToHT(mdl, state)
%INPUTS
%mdl:           WmrModel object
%state:         ns x 1 state vector
%OUTPUTS
%HT_to_parent:  4x4 x nf Homogeneous transform from each frame to parent frame
%HT_to_world:   (optional) 4x4 x nf Homogeneous transform from each frame to world frame


nf = mdl.nf;
ns = length(state);
[isorient,ispos,isjd] = stateIs(ns);

%TODO, eliminate this concatenation?
HT_parent_jd0 = cat(3,mdl.frames.HT_parent_jd0); 
dof_type = [mdl.frames.dof_type];


HT_parent = HT_parent_jd0; %init to _jd0 so that
%only R must be updated for revolute joints, 
%only t must be update for prismatic joints
HT_parent(:,:,1) = poseToHT(state(isorient),state(ispos));

%COMPUTE HT_PARENT

%this function should be called whenever '.joint_pos' changes
%I is a list of indices i;

jds = state(isjd); %joint displacements

for fi = 2:nf
    
    jd = jds(fi-1); %joint displacement

    %HT^(parent)_(child) = HT^(parent)_(child,jp=0)*HT^(child,jp=0)_(child)

    %for revolute joint about x:
    % = [R, t; 0 0 0 1] * [Rx(jp), [0 0 0]'; 0 0 0 1]; 
    % = [R*Rx(jp), t; 0 0 0 1];
    %rotate about the x axis of frame i, not the x axis of the parent frame

    %for prismatic joint along x:
    % = [R, t; 0 0 0 1] * [eye(3), jp*[1 0 0]'; 0 0 0 1]; 
    % = [R, R*(jp*[1 0 0]')+t; 0 0 0 1];
    %translate along the x axis of frame i, not the x axis of the parent frame

    dof_type_ = dof_type(fi);
    if dof_type_ < 4
        %revolute
        axis_no = dof_type_;
        if axis_no == 1
            Rota = Rotx(jd);
        elseif axis_no == 2
            Rota = Roty(jd);
        elseif axis_no == 3
            Rota = Rotz(jd);
        end
        
        R = HT_parent_jd0(1:3,1:3,fi);
        HT_parent(1:3,1:3,fi) = R*Rota;
    else
        %prismatic
        axis_no = dof_type_ - 3;
        r = HT_parent_jd0(axis_no,1:3,fi);
        t = HT_parent_jd0(1:3,4,fi);
        HT_parent(1:3,4,fi) = t+r'*jd;
    end

end


%COMPUTE HT_WORLD
if nargout > 1
    %HT^(world)_(frame i) = HT^(world)_(parent)*HT^(parent)_(frame i)
    parent_ind = [mdl.frames.parent_ind];

    HT_world = zeros(4,4,mdl.nf);
    HT_world(:,:,1) = HT_parent(:,:,1);
    for fi = 2:nf
        HT_world(:,:,fi) = HT_world(:,:,parent_ind(fi)) * HT_parent(:,:,fi);
    end
end




