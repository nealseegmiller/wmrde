function contacts = updateWheelContactGeom(surfaces, HT_wheel_to_world, radii)
%Updates wheel-ground contact geometry.
%Discretizes the surface of each wheel into points, and chooses the point for which dz is smallest
%INPUTS
%surfaces:              cell array of surface objects
%HT_wheel_to_world:     4x4 x nw,  Homogeneous transforms from wheel to world coords
%radii:                 1 x nw, wheel radii
%OUTPUT
%contacts               1 x nw array of WheelContactGeom objects


nw = length(radii);


%no need check for contact on the top half of the wheel

%projection of world z axis onto wheel x axis. 
%px(i) = [0 0 1]*HT_wheel_to_world(1:3,1,i)
px = zeros(1,nw); 
px(:) = HT_wheel_to_world(3,1,:);

%projection of world z axis onto wheel z axis. 
%pz(i) = [0 0 1]*HT_wheel_to_world(1:3,3,i)
pz = zeros(1,nw); 
pz(:) = HT_wheel_to_world(3,3,:);

mid = atan2(px,pz); %middle angle


rng = pi;
np = 18+1;
% np = 180*10+1; %DEBUGGING

%potential contact points for all wheels
%evenly distributed along arc of wheel circumference
[allcp,allcp_world,allangles] = discretizeToPoints(radii,mid,rng,np,HT_wheel_to_world);


if 1 %refine contact point location with 2nd step (smaller range, finer discretization)
    %TODO, make the 2nd step faster by using the surface index from the 1st step?

    rng = rng/(np-1);
    np = 10+1;

    [~,~,iall] = minDz(surfaces,allcp_world,nw);
    mid = allangles(iall);
    [allcp,allcp_world,allangles] = discretizeToPoints(radii,mid,rng,np,HT_wheel_to_world);
end     

[dz,si,iall] = minDz(surfaces,allcp_world,nw);

cp = allcp(:,iall);
cp_world = allcp_world(:,iall);
contact_angles = allangles(iall);
N_world = surfacesNormal(surfaces,cp_world,si);


%copy to WheelContactGeom
contacts(1:nw) = WheelContactGeom();
for wno = 1:nw
    contacts(wno).dz = dz(wno);
    contacts(wno).angle = contact_angles(wno);
    N = HT_wheel_to_world(1:3,1:3,wno)' * N_world(:,wno); %in wheel coords
    contacts(wno).HT_wheel = HTContactToWheel(cp(:,wno),N);
    contacts(wno).HT_world = HT_wheel_to_world(:,:,wno) * contacts(wno).HT_wheel;
end

end


function [allcp, allcp_world, allangles] = discretizeToPoints(radii,mid,rng,np,HT_wheel_to_world)
    %given an arc of the wheel circumference specified by mid,rng
    %discretize to possible contact points
    %cat points for *all* wheels into a single matrix to minimize calls to surfaces functions
    
    nw = length(radii);
    
    angles = cell(1,nw);
    cp = cell(1,nw);
    cp_world = cell(1,nw);
    
    da = rng/(np-1);
    for wno = 1:nw
        angles{wno} = (mid(wno)-rng/2):da:(mid(wno)+rng/2);
        cp{wno} = contactAnglesToPoints(radii(wno),angles{wno});
        cp_world{wno} = applyHT(HT_wheel_to_world(:,:,wno),cp{wno});
    end
    
    allangles = [angles{:}];
    allcp = [cp{:}];
    allcp_world = [cp_world{:}];
    
end


function [dz,si,iall] = minDz(surfaces,allcp_world,nw)
    %choose the contact point for each wheel with minimum dz
    
    dz = zeros(1,nw); %contact height errors
    si = zeros(1,nw); %surface indices
    iall = zeros(1,nw); %indices in allcp
    
    np = size(allcp_world,2)/nw; %number of points per wheel
    whichw = uncatIndexArray(np*ones(1,nw)); %which wheel contact point belongs to
    
    [alldz,allsi] = surfacesDz(surfaces,allcp_world); %just one call
    for wno = 1:nw
        [~,i] = min(alldz(whichw == wno));
        iall_ = i + (wno-1)*np;
        
        dz(wno) = alldz(iall_);
        si(wno) = allsi(iall_);
        iall(wno) = iall_;
    end

end




