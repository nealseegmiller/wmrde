function contacts = updateWheelContactGeomBlind(HT_wheel_to_world, radii, contact_angles)
%updates wheel contact geometry without perception of terrain surface. Requires contact angles as input.
%INPUTS
%HT_wheel_to_world: 4x4 x nw, Homogeneous transform from wheel to world coords
%radii:             1 x nw, wheel radii
%contact_angles:    1 x nw
%OUTPUT
%contacts:          1 x nw array of WheelContactGeom objects

nw = length(radii);

cp = zeros(3,nw); %contact points in wheel coords
N = zeros(3,nw); %surface normals in wheel coords

for wno = 1:nw
    cp(:,wno) = contactAnglesToPoints(radii(wno),contact_angles(wno));
    %assume surface normal points from contact point to wheel frame origin
    N_ = -cp(:,wno);
    N_ = N_/norm(N_); %normalize
    N(:,wno) = N_;
end

%copy to WheelContactGeom
contacts(1:nw) = WheelContactGeom();
for wno = 1:nw
    contacts(wno).angle = contact_angles(wno);
    contacts(wno).HT_wheel = HTContactToWheel(cp(:,wno),N(:,wno));
    contacts(wno).HT_world = HT_wheel_to_world(:,:,wno) * contacts(wno).HT_wheel;
end


