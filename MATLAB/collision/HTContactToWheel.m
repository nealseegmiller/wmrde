function HT_contact_to_wheel = HTContactToWheel(cp,N)
%INPUTS
%cp:    3 x nw, contact points (in wheel coords)
%N:     3 x nw, terrain normal vectors (in wheel coords)

nw = size(cp,2);
HT_contact_to_wheel = zeros(4,4,nw);
HT_contact_to_wheel(4,4,:) = 1;

for wno = 1:nw

    zaxis = N(:,wno); %z axis = surface normal

    xaxis = cross3([0 1 0]',zaxis); %x axis = cross(wheel frame y axis, z axis)
    xaxis = xaxis/norm(xaxis); %normalize to unit vector
    yaxis = cross3(zaxis,xaxis);
    
    HT_contact_to_wheel(1:3,:,wno) = [xaxis, yaxis, zaxis, cp(:,wno)];
end