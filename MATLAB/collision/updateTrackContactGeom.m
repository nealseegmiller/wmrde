function contacts = updateTrackContactGeom(surfaces, HT_track_to_world, contacts)

nt = length(contacts); %number of tracks

cp_world = cell(1,nt); %contact points in world coords
Nz_world = cell(1,nt); %z components of normal vectors in world coords

for tno = 1:nt
    np = size(contacts(tno).HT_track,3);
    
    %HT_contact_to_world
    contacts(tno).HT_world = zeros(4,4,np);
    for pno = 1:np
        contacts(tno).HT_world(:,:,pno) = HT_track_to_world(:,:,tno) * contacts(tno).HT_track(:,:,pno);
    end
    
    cp_world{tno} = zeros(3,np);
    cp_world{tno}(:) = contacts(tno).HT_world(1:3,4,:);
    
    tmp = contacts(tno).HT_world(3,3,:);
    Nz_world{tno} = tmp(:)';
    
end

allcp_world = [cp_world{:}];
allup = [Nz_world{:}] > 0; %true if normal points up

allzsurf = zeros(1,length(allup));
allzsurf(allup) = surfacesHeight(surfaces,allcp_world(:,allup)); %single call

whichtrack = uncatIndexArray([contacts.np]); %which track each point belongs to

for tno = 1:nt
    zsurf = allzsurf(whichtrack == tno);
    dh = cp_world{tno}(3,:) - zsurf; %delta height
    
    up = Nz_world{tno} > 0;
    contacts(tno).dz(up) = dh(up) .* Nz_world{tno}(up);
    contacts(tno).dz(~up) = Inf; %don't use these
end




