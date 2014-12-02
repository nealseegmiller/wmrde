% clear_all
% close all

% %initialize random number stream
% %RandStream.list for available gentypes
% stream = RandStream.create('mt19937ar','seed',147);
% RandStream.setGlobalStream(stream);


function surfaces = rocks(surfaces)
%makes rock surface objects
%TODO, 
%-detail rocks
%-make triangle size more uniform for oblong rocks?

%OPTIONS
nrocks = 10;

nv = 20;
roundness = .1;

%rock size
lx_min = 2;
lx_max = 4;
ly_min = 2;
ly_max = 4;

dsplit = 1;
npert = .1;

%center of boulder field
cx = 4;
cy = 0;
%range of boulder field
rx = 2;
ry = 9;
%range of yaw values
ryaw = pi;

%END OPTIONS

L = ones(3,nrocks);
L(1,:) = (rand(1,nrocks)-.5)*(lx_max-lx_min) + (lx_max + lx_min)/2;
L(2,:) = (rand(1,nrocks)-.5)*(ly_max-ly_min) + (ly_max + ly_min)/2;
L(3,:) = min(L(1:2,:));

L(:,3:end) = L(:,3:end)/3; %DEBUGGING

%get rock vertices
V = cell(1,nrocks);
I = cell(1,nrocks);
for i = 1:nrocks
    [V{i},I{i}] = rock(nv,roundness,L(:,i),dsplit,npert);
end

%rock orientation, position
yaw = zeros(1,nrocks);
C = zeros(3,nrocks); %rock center locations
VR = cell(1,nrocks); %vertices of bounding rectangle for each rock
nr = 0;

while 1
    
    %get vertices of bounding rectangle
    V_ = V{nr+1};
    
    xmin = min(V_(1,:));
    xmax = max(V_(1,:));
    ymin = min(V_(2,:));
    ymax = max(V_(2,:));
    
    VR_ = [xmin xmin xmax xmax;
           ymin ymax ymax ymin;
           zeros(1,4)];

    %rotate
    yaw_ = (rand-.5)*ryaw;
    VR_ = Rotz(yaw_)*VR_; 
    
    %translate
    x = (rand-.5)*rx + cx;
    y = (rand-.5)*ry + cy;
    VR_(1,:) = VR_(1,:) + x;
    VR_(2,:) = VR_(2,:) + y;
    
    %close the polygon
    VR_(:,5) = VR_(:,1);
    
    collision = false;
    for i = 1:nr
        [xi,yi] = polybool('intersection', VR_(1,:),VR_(2,:),VR{i}(1,:),VR{i}(2,:)); %slow?
        if ~isempty(xi)
            collision = true;
            break
        end
    end
    
    if collision == false
        nr = nr+1;
        yaw(nr) = yaw_;
        C(:,nr) = [x y 0]';
        VR{nr} = VR_;
    end
    
    if nr >= nrocks
        break
    end
    
end

%bury rocks halfway into ground
C(3,:) = surfacesHeight(surfaces,C(1:2,:));


%make TriMeshSurf objects
surfaces = cell(1,nrocks);
for i = 1:nrocks
    V_ = V{i};
    I_ = I{i};
    
    %rotate
    Rz = Rotz(yaw(i));
    V_ = Rz*V_;
    
    %translate
    nv_ = size(V_,2);
    V_ = V_ + C(:,i)*ones(1,nv_);
    
    surfaces{i} = TriMeshSurf();
    
    if 0
        addVertices(surfaces{i},V_);
        addTriangles(surfaces{i},I_);
    
    else
        %duplicate vertices to fix lighting
        nt = size(I_,1);
        IT = I_';

        addVertices(surfaces{i},V_(:,IT(:)));
        addTriangles(surfaces{i},reshape(1:3*nt,3,nt)');
    end

end


if 1
    %%
    %PLOT
    
    set(figure,'name','rocks')
    hold on
    drawSurfaces(surfaces);
    
    axis equal
    
    xlabel('x')
    ylabel('y')
    zlabel('z')
    makeLegible(14)
    
%     h_light = light('Position',[1 0 1],'Parent',gca);
    
    %DEBUGGING
    for i = 1:nrocks
        plot(VR{i}(1,:),VR{i}(2,:))
    end
end

