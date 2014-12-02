function pts = contactPointsOnSurface(contacts)
%returns the wheel-terrain contact point on the terrain surface
%INPUTS
%contacts:  1xnw array of WheelContactGeom objects

nw = length(contacts);
pts = NaN(3,nw);

for i = 1:nw
    if contacts(i).dz < 0
        normal = contacts(i).HT_world(1:3,3);
        pts(:,i) = contacts(i).HT_world(1:3,4) - normal*contacts(i).dz;
    end
end
