function [meas,inpt] = readdataLandtamer(filenames,starttime,stoptime,p,cinfo,mdl)


opts.downsample=1;
opts.plot=1;

if length(filenames) > 1
    
    i = strfindCell(filenames,'_commands.data');
    commands_load = load(filenames{i});
    %Time stamp, Actuation Time, Velocity, Curvature
    
    i = strfindCell(filenames,'_pose.data');
    pose_load = load(filenames{i});
    %Time stamp, X, Y, Z, Roll, Pitch, Yaw, Velocity, Curvature

    i = strfindCell(filenames,'_state.data');
    state_load = load(filenames{i});
    % System Time, time_sec, time_usec, state, enc_left, enc_right, tca_left,
    % tca_right, tach, req_lvel, req_rvel, cur_lvel, cur_rvel, p, i, d, acc,
    % throttle
    
    %save to .mat file for faster loading
    matfilename=strrep(filenames{1},'_commands.data','.mat');
    save(matfilename,'commands_load','pose_load','state_load');
    
else
    %load from .mat file, faster!
    load(filenames{1})
end

pose = pose_load;
commands = commands_load;
state = state_load;

%meas
meas.t = pose_load(:,1);

% blah=1;

if 0
    %DEBUGGING, FOR FSR 2013
    nm = length(meas.t);
    %train
%     starttime = meas.t(round(nm/4));
%     stoptime = meas.t(round(nm/2));
    %test
    starttime = meas.t(3400);
    stoptime = meas.t(round(nm/4));
end

if ~isempty(starttime) 

    startind = find(meas.t > starttime,1);
    stopind = find(meas.t < stoptime,1,'last');
    
    pose = pose(startind:stopind,:);
    meas.t = meas.t(startind:stopind);
end


if opts.downsample
    keepevery = 4;
    meas.t = meas.t(1:keepevery:end);
    pose = pose(1:keepevery:end,:);
else

    %remove rows if no gps update
    %not necessary if downsample
    noupdate = all(numdiff(pose(:,2:4),[],-1)==0,2);
    meas.t = meas.t(~noupdate);
    pose = pose(~noupdate,:);
end

nm=length(meas.t);

meas.pos = pose(:,2:4);

%TODO, eliminate rows where no gps update, meas.pos == circshift(meas.pos,1)
%causes velocity error

%handle bad gps
d = sqrt(sum((ones(nm,1)*meas.pos(1,:) - meas.pos).^2,2)); %distance from starting position
isbad=d>1e4;
meas.pos(isbad,:) = interp1(meas.t(~isbad), meas.pos(~isbad,:), meas.t(isbad), 'linear');

%start at [0 0 0]. numerical error in animation otherwise
gps_start = meas.pos(1,:);

if 1
    %DEBUGGING
    hemisphere = 'N';
    lonzone = 17; %Pittsburgh: 17,

    [meas.lat,meas.lon] = UTMToLatLon(meas.pos(:,2),meas.pos(:,1),hemisphere,lonzone);
end


meas.pos = meas.pos - ones(nm,1)*gps_start; 



%distance traveled
meas.s = zeros(nm,1);
for i = 2:nm
    meas.s(i) = meas.s(i-1) + norm(meas.pos(i,:) - meas.pos(i-1,:));
end

meas.euler = pose(:,5:7);
meas.euler(:,3) = unwrapRad(meas.euler(:,3)); %unwrap yaw




%inpt
inpt.t = state(:,1); %???
% inpt.t = state(:,2) + state(:,3)*1e-6;

if ~isempty(starttime)
    startind = find(inpt.t > starttime,1);
    stopind = find(inpt.t < stoptime,1,'last');
    
    state = state(startind:stopind,:);
    inpt.t = inpt.t(startind:stopind);
    
end

%start at time = 0
t0 = min(meas.t(1),inpt.t(1));
meas.t = meas.t - t0;
inpt.t = inpt.t - t0;


% intom = 2.54/100; %inch to meter
% rad=33/2*intom; %wheel radius

radii = [mdl.frames(mdl.wheelframeinds).rad];
rad = radii(1);

inpt.enc = state(:,[12 13 12 13 12 13])/rad;

inpt.rol=interp1(meas.t,meas.euler(:,1),inpt.t,'nearest','extrap');
inpt.pit=interp1(meas.t,meas.euler(:,2),inpt.t,'nearest','extrap');



if opts.plot    
    
    %%
    
    if 1
        %%
        zoom = 20;
        border = 2;
        [h_im,h_line] = plotPathOnAerial(meas.pos(:,1),meas.pos(:,2),meas.lat,meas.lon,zoom,border);
        plot(meas.pos(1,1),meas.pos(1,2),'b>')
        plot(meas.pos(end,1),meas.pos(end,2),'bs')
    end
    
    set(figure,'name','path 3D')
    hold on
    plot3(meas.pos(:,1), meas.pos(:,2), meas.pos(:,3),'DisplayName','meas')
    xlabel('x')
    ylabel('y')
    zlabel('z')
    legend('Location','Best')
    axis equal
    makeLegible(14)
    
    set(figure,'name','orientation')
    plot(meas.t,meas.euler*180/pi)
    legend('rol','pit','yaw')
    xlabel('time (s)')
    ylabel('deg')
    makeLegible(14)
    
end


end








