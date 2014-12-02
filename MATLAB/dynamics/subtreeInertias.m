function Is_subt = subtreeInertias(mdl,Xup)
%INPUT
%mdl:   WmrModel object
%Xup:   6x6xnf Plucker transforms, Xup(:,:,i)' transforms spatial *force* vector from frame i to parent(i) coords
%OUTPUT
%Is:    6x6xnf, Is(:,:,i) is the 6x6 spatial inertia of the subtree rooted at frame i, in the coordinates of frame i

parent_ind = [mdl.frames.parent_ind];

%initialize to root frame inertia
Is_subt = cat(3,mdl.frames.Is);
for i = mdl.nf:-1:2
    %transform & add inertia to parent frame
    %even if a parent frame has no mass, the subtree rooted at that frame may already have mass due to other children
    j = parent_ind(i);
    Is_subt(:,:,j) = Is_subt(:,:,j) + Xup(:,:,i)'*Is_subt(:,:,i)*Xup(:,:,i);
end