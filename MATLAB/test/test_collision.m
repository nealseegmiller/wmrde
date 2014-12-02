%test_collision.m
% test Surface class and collision detection functions

%%
%test_surface
clear_all
close all
clc

surfs = {};

% surfs{end+1} = flat;
% surfs{end+1} = ramp();

setseed(123);
surfs{end+1} = randomgrid();
printGridSurf(surfs{end},[resourcedir() 'gridsurfdata.txt'])


lo = 2.1;
hi = 4.1;
np = 5;

d = (hi-lo)/max([np-1, 1]);

pts = zeros(3,np);
for i = 1:np
    pts(1,i) = lo + d*(i-1);
    pts(2,i) = -1;
    pts(3,i) = .25;
end


pts

height = surfacesHeight(surfs,pts(1:2,:))
[dz,surface_index] = surfacesDz(surfs,pts)
normal = surfacesNormal(surfs,pts(1:2,:),surface_index)


if 1
    %time it
    irf = 1e-2;
    n = 1e5;    
    
    tic
    for i=1:(n*irf)
        [dz,~] = surfacesDz(surfs,pts);
    end
    t=toc;
    
    fprintf('number of points: %d\n',np)
    fprintf('iterations: %e\n',n)
    t_=t/irf*1e3;
    fprintf('Elapsed time (ms): %f\n',t_)

    
end

if 1
    %%
    %PLOT
    set(figure, 'name','surf')
    hold on
    drawSurfaces(surfs,gca);

    plot3(pts(1,:),pts(2,:),pts(3,:),'b.')
    
    uvw = normal .* -(ones(3,1)*dz);
    h = drawQuivers(pts,uvw);
    axis equal
    
    d = .5;
    lim = [min(pts,[],2) max(pts,[],2)]';
    lim = lim(:)' + [-d d -d d -d/2 d/2];
    axis(lim)
    
    view(30,30);
    
    xlabel('x')
    ylabel('y')
    zlabel('z')
    makeLegible(14)
    
    resizeFig(1.5,1.5)
end



%%
%test_updateWheelContactGeom()

clear_all
close all
clc

rad=.325;

%terrain
surfs = {};

% surfs{end+1} = flat();
% surfs{end+1} = ramp();

setseed(321);
surfs{end+1} = randomgrid();
printGridSurf(surfs{end},[resourcedir() 'gridsurfdata.txt'])

%wheel pose
orient = [0 0 0]'*pi/180;
pos = [2.5,-1,rad]';
pos(3) = pos(3) + surfacesHeight(surfs,pos(1:2));

HT_wheel_to_world = poseToHT(orient,pos);

contact = updateWheelContactGeom(surfs, HT_wheel_to_world, rad);

disp('wheel pos = ')
disp(pos)
disp(['dz = ' num2str(contact.dz)])
disp(['contact_angle (deg) = ' num2str(contact.angle*180/pi)])
disp('HT_contact_to_wheel = ')
disp(contact.HT_wheel)
disp('HT_contact_to_world = ')
disp(contact.HT_world)



if 1
        
    %time it
    irf=1e-3;
    n=1e5;

    %vectorize
    nw=4; %number of wheels
    radii=rad*ones(1,nw);
    HT_wheels_to_world = repmat(HT_wheel_to_world,[1 1 nw]);

    % profile on
    tic
    for i=1:(n*irf);
        blah = updateWheelContactGeom(surfs, HT_wheels_to_world, radii);
    end
    t=toc;
    % profile off
    % profile viewer

    fprintf('num wheels: %d\n',nw)
    fprintf('iterations: %e\n',n)
    t_=t/irf*1e3;
    fprintf('Elapsed time (ms): %f\n',t_)
end

if 1
    %%
    %PLOT
    set(figure, 'name','surf')
    hold on
    drawSurfaces(surfs,gca);
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('z')
    makeLegible(14)
    
    %wheel
    anim = hgtransform('Parent',gca);
    nv = 36;
    C = 'none';
    makeCircleHgt(anim,rad,2,zeros(3,1),nv,C);
    set(anim,'Matrix',HT_wheel_to_world)
    
    %normal
    hq = drawQuivers(contact.HT_world(1:3,4),contact.HT_world(1:3,3)*rad*.75);
    set(hq,'Color','g','MaxHeadSize',.5,'LineWidth',1.5)
    
    %center axes on wheel
    c = HT_wheel_to_world(1:3,4)';
    c = [c(1) c(1) c(2) c(2) c(3) c(3)];
    d = 5*rad;
    axislimits = [-d d -d d -d/2 d/2];
    
    axis(c+axislimits)
    
    view(0,0);
%     resize_fig(2,2)

    drawnow

    
end


if 1
    %%
    %DEBUGGING, validate ddz/dca calculation
%     close all
    
    np = 180*10+1;
    angles = linspace(-pi/2,pi/2,np);
    pts_wheel = contactAnglesToPoints(rad,angles);
    pts_world = applyHT(HT_wheel_to_world,pts_wheel);

    [dz,si] = surfacesDz(surfs,pts_world);
    N = surfacesNormal(surfs,pts_world,si);

    
    %tangent - cos(a)*x axis + sin(a)*z axis
    T = zeros(3,np);
    for i = 1:np
        T(:,i) = - cos(angles(i))*HT_wheel_to_world(1:3,1) + sin(angles(i))*HT_wheel_to_world(1:3,3);
    end
    
    ddz_dca_pred = sum(T.*N,1)*rad;  %predicted

    ddz_dca = numDiff(dz',angles',0)'; %actual
    
    
    %PLOT
    
    set(figure,'name','ddz/dca')
    hold on
    plot(angles*180/pi,ddz_dca,'b:')
    plot(angles*180/pi,ddz_dca_pred,'g')
    plot([angles(1) angles(end)]*180/pi,[0 0],'r')
    myaxis=axis;
    xlabel('contact angle (deg)')
    ylabel('m/rad')
    legend('act','pred','Location','Best')
    makeLegible(14)
    
    
    set(figure, 'name','dz')
    hold on
    plot(angles*180/pi,dz,'b')
    clear h
    [~,i] = min(dz);
    h(1) = plot(angles(i)*180/pi,dz(i),'bo','DisplayName',sprintf('act: %.1f deg',angles(i)*180/pi));
    [~,i] = min(abs(ddz_dca_pred)); %TODO, better rootfinding
    h(2) = plot(angles(i)*180/pi,dz(i),'go','DisplayName',sprintf('pred: %.1f deg',angles(i)*180/pi));
    
    
    xlabel('contact angle (deg)')
    ylabel('dz (m)')
    legend(h,'Location','Best')
    
    makeLegible(14)


    
    
end

return

%%
%test_updateTrackContactGeom()

clear_all
close all
clc

%make track frame
frame = Frame();
frame.rad = .2;
frame.rad2 = .2;
frame.L = .5;

contact = initTrackContactGeom(frame);

%terrain
surfs = {};
% surfs{end+1} = flat();
% surfs{end+1} = ramp();

setseed(123);
surfs{end+1} = randomgrid();
printGridSurf(surfs{end},[resourcedir() 'gridsurfdata.txt'])

%track pose
orient = [0 0 0]'*pi/180;
pos = [3 0 frame.rad - .05]';
pos(3) = pos(3) + surfacesHeight(surfs,pos)
HT_track_to_world = poseToHT(orient,pos);

contact = updateTrackContactGeom(surfs, HT_track_to_world, contact);

contact.dz
contact.dvc_dvt

if 1
    %%
    %PLOT
    set(figure, 'name','surf')
    hold on
    drawSurfaces(surfs,gca);
    
    np = size(contact.HT_world,3);
    
    %plot contact points
    cp = zeros(3,np);
    cp(:) = contact.HT_world(1:3,4,:);
    
    plot3(cp(1,:),cp(2,:),cp(3,:),'b.-')
    
    %quivers from contact points to surface
    N = zeros(3,np);
    N(:) = contact.HT_world(1:3,3,:);
    uvw = -N .* (ones(3,1)*contact.dz);
    
    incontact = contact.dz < 0;

    h = drawQuivers(cp(:,incontact),uvw(:,incontact));
    set(h,'Color','green')

    axis equal
    
    d = .5;
    lim = [min(cp(:,incontact),[],2) max(cp(:,incontact),[],2)]';
    lim = lim(:)' + [-d d -d d -d d];
    axis(lim);
    
    view(30,30)
    
    xlabel('x')
    ylabel('y')
    zlabel('z')
    makeLegible(14)
    
    resizeFig(1.5,1.5)
    
end
