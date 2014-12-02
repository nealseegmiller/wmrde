function [HT_contact_to_track] = HTContactToTrack(rad1,rad2,L,npflat,sides)
%Get equally spaced points along track surface, in track frame coordinates
%track frame:
%-origin is located at center of 1st sprocket
%-x axis direction is toward center of 2nd sprocket
%-y axis is sprocket rotation axis
%INPUTS
%rad1:      radius of 1st sprocket
%rad2:      radius of 2nd sprocket
%L:         distance between sprockets (may be negative)
%npflat:    number of points on flat section
%sides:     1x4 logical, [lower flat, outer sprocket, upper flat, inner sprocket]
%OUTPUT
%HT_contact_to_track: 4 x 4 x np, homogeneous transform from contact to track coords for each point
%           contact frame x axes are tangent to terrain surface, y axes are parallel to track frame y axis


HT_contact_to_track = zeros(4,4,0);

ca = -asin((rad1-rad2)/L); %contact angle
Lflat = L*cos(ca);
dx = Lflat/(npflat-1);

if sides(1) %lower flat
    
    %first contact point
    cp = contactAnglesToPoints(rad1,ca); 
    N = -cp;
    N = N/norm(N);
    HT = HTContactToWheel(cp,N);

    %all contact points
    cp = HT(1:3,4)*ones(1,npflat) + HT(1:3,1)*(0:dx:Lflat);
    HT = HTContactToWheel(cp,N*ones(1,npflat));
    
    if sides(4) %avoid duplicate points
        HT(:,:,1) = [];
    end
    
    HT_contact_to_track = cat(3,HT_contact_to_track,HT);
end

rca = sign(L)*pi + 2*ca; %range of contact angles for front sprocket

if sides(2) %front sprocket
    npsprocket = round(rca*rad2/dx) + 1;
    dca = rca/(npsprocket-1);
    angles = ca:-dca:(ca-rca);
    
    cp = contactAnglesToPoints(rad2,angles);
    cp(1,:) = cp(1,:) + L;
    
    N = zeros(3,npsprocket);
    for pno = 1:npsprocket
        N_ = [L 0 0]' - cp(:,pno);
        N_ = N_/norm(N_);
        N(:,pno) = N_;
    end
    
    HT = HTContactToWheel(cp,N);
    
    if sides(1) %avoid duplicate points
        HT(:,:,1) = [];
    end
    
    HT_contact_to_track = cat(3,HT_contact_to_track,HT);
end


if sides(3) %upper flat
    %first contact point
    cp = contactAnglesToPoints(rad2,ca-rca); 
    cp(1,:) = cp(1,:) + L;
    N = [L 0 0]' - cp;
    N = N/norm(N);
    HT = HTContactToWheel(cp,N);

    %all contact points
    cp = HT(1:3,4)*ones(1,npflat) + HT(1:3,1)*(0:dx:Lflat);
    HT = HTContactToWheel(cp,N*ones(1,npflat));
    
    if sides(2) %avoid duplicate points
        HT(:,:,1) = [];
    end
    
    HT_contact_to_track = cat(3,HT_contact_to_track,HT);    
end


rca = sign(L)*pi - 2*ca; %range of contact angles for rear sprocket
if sides(4) %rear sprocket
    npsprocket = round(rca*rad1/dx) + 1;
    dca = rca/(npsprocket-1);
    angles = (ca+rca):-dca:ca;
    
    cp = contactAnglesToPoints(rad1,angles);
    
    N = zeros(3,npsprocket);
    for pno = 1:npsprocket
        N_ = -cp(:,pno);
        N_ = N_/norm(N_);
        N(:,pno) = N_;
    end
    
    HT = HTContactToWheel(cp,N);
    
    if sides(3) %avoid duplicate points
        HT(:,:,1) = [];
    end
    
    HT_contact_to_track = cat(3,HT_contact_to_track,HT);
end





