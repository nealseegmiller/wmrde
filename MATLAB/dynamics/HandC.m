function [H, inv_H, C] = HandC(mdl, HT_parent, qvel)
%Compute the joint space inertia and bias force for forward dynamics

%OPTIONS
do_reuse_H = 0;
do_approx_C = 0; %treat WMR as single rigid body when computing joint space bias force

%TODO, eliminate persistent vars, store in WmrModel instead
persistent H_prev inv_H_prev Ist

%spatial transforms
Xup = invHTsToPluckers(HT_parent);
%Xup(:,:,i) transforms spatial *motion* vector from frame parent(i) to i coords
%Xup(:,:,i)' transforms spatial *force* vector from frame i to parent(i) coords!

%joint space inertia
if ~do_reuse_H || isempty(H_prev)

    %Ist(:,:,i) is the 6x6 spatial inertia of the subtree rooted at frame i, in the coordinates of frame i
    Ist = subtreeInertias(mdl,Xup);
    H = jointSpaceInertia(mdl,Xup,Ist);
    
    H_prev = H; 
else
    H = H_prev;
end

if mdl.use_constraints
    if ~do_reuse_H || isempty(inv_H_prev)
        %invert H
        %matrix inversion is slow & innaccurate in general, but inv(H) is required many times so precomputing saves time?
        %faster than inv(H)?

        U = chol(H);
        invU = U\eye(size(U)); %backslash operator recognizes upper triangular systems
        inv_H = invU*invU';

        inv_H_prev = inv_H;
    else
        inv_H = inv_H_prev;
    end
else
    inv_H = []; %not necessary
end

%joint space bias force
if ~do_approx_C
    C = jointSpaceBiasForce(mdl,Xup,qvel);
else
    %approximate, all inertia at body frame
    %spatial acc of gravity in body coords
    g = zeros(6,1);
    g(4:6) = Xup(1:3,3,1)*-mdl.grav; %world z axis in body coords * grav acc

    C = zeros(nv,1);
    C(1:6) = -Ist(:,:,1)*g + crf(qvel(1:6))*Ist(:,:,1)*qvel(1:6);
end


