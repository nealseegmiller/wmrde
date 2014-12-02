function J = DvcontactDvtrack(HT_contact_to_track, rad)
%HT_contact_to_track:   4x4 x np, homogeneous transform from contact to track coords for each point
%rad:                   radius of driving sprocket for track
%J:                     3*np x 7, Jacobian of contact point linear vel wrt [track frame spatial vel; sprocket angular rate]

%Jacobian of contact point linear vel ...
np = size(HT_contact_to_track,3);
J = zeros(3*np,7);

%contact point (in track coords)
cp = zeros(3,np);
cp(:) = HT_contact_to_track(1:3,4,:);
%contact frame x axes (in track coords)
X = zeros(3,np);
X(:) = HT_contact_to_track(1:3,1,:);

%for indexing spatial vectors
ang = 1:3;
lin = 4:6;

for pno = 1:np

    I = (pno-1)*3 + (1:3); %row indices

    %... wrt track frame spatial vel
    J(I,ang) = skew3(cp(:,pno))';
    J(I,lin) = eye(3);

    %... wrt sprocket angular rate
    J(I,7) = -X(:,pno)*rad;

    %rotate into contact frame coords
    J(I,:) = HT_contact_to_track(1:3,1:3,pno)' * J(I,:);

end