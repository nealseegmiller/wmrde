function [meas,inpt] = readdataRecbot(filenames,starttime,stoptime,p,cinfo,mdl)

opts.downsample = 1;
opts.plot = 1;
opts.smooth = 1;
opts.inityaw = 1;

if length(filenames) > 1
    novatel_load = load(filenames{1});
    %Time stamp, X, Y, Z, Roll, Pitch, Yaw
    
    status_load = load(filenames{2});
    %Time stamp, Commanded speed, Turn angle encoder
    
    %save to .mat file for faster loading
    matfilename=strrep(filenames{1},'_novatel.data','.mat');
    save(matfilename,'novatel_load','status_load');
    
else
    %load from .mat file, faster!
    load(filenames{1})
end

pose=novatel_load;
state=status_load;

%meas
meas.t = pose(:,1);

% blah=1;

if ~isempty(starttime) 

    startind = find(meas.t > starttime,1);
    stopind = find(meas.t < stoptime,1,'last');
    
    pose = pose(startind:stopind,:);
    meas.t = meas.t(startind:stopind);
end


if opts.downsample
    keepevery=10;
    meas.t = meas.t(1:keepevery:end);
    pose = pose(1:keepevery:end,:);
else

    %remove rows if no gps update
    %not necessary if downsample
    pose_cs = circshift(pose,1);
    noupdate = all(pose(:,2:4)==pose_cs(:,2:4),2);
    meas.t = meas.t(~noupdate);
    pose = pose(~noupdate,:);
end

nm=length(meas.t);

meas.pos = pose(:,2:4);

%TODO, eliminate rows where no gps update, meas.pos == circshift(meas.pos,1)
%causes velocity error

%handle bad gps
% d = sqrt(sum((ones(nm,1)*meas.pos(1,:) - meas.pos).^2,2)); %Euclidean distance from starting position
% isbad=d>1e4;
% meas.pos(isbad,:) = interp1(meas.t(~isbad), meas.pos(~isbad,:), meas.t(isbad), 'linear');

%start at [0 0 0]. numerical error in animation otherwise
gps_start = meas.pos(1,:);

if 1
    %DEBUGGING
    hemisphere = 'N';
    lonzone = 17; %Pittsburgh: 17
    
    [meas.lat,meas.lon] = UTMToLatLon(meas.pos(:,2),meas.pos(:,1),hemisphere,lonzone);

end


meas.pos = meas.pos - ones(nm,1)*gps_start; 

%distance traveled
meas.s = zeros(nm,1);
for i = 2:nm
    meas.s(i) = meas.s(i-1) + norm(meas.pos(i,:) - meas.pos(i-1,:));
end

meas.euler = pose(:,5:7)*pi/180; %degrees to radians
meas.euler(:,3) = unwrapRad(meas.euler(:,3)); %unwrap yaw


%inpt
inpt.t = state(:,1);

if ~isempty(starttime)
    startind = find(inpt.t > starttime,1);
    stopind = find(inpt.t < stoptime,1,'last');
    
    state = state(startind:stopind,:);
    inpt.t = inpt.t(startind:stopind);
    
end

if opts.downsample
    keepevery=4;
    inpt.t = inpt.t(1:keepevery:end);
    state = state(1:keepevery:end,:);
end

meas.omit = false(nm,1);
if 1
    bad = meas.t >= 417235.3 & meas.t <= 417237.1;
    meas.omit(bad) = true;
    if any(bad)
        disp([mfilename ': omit due to bad speed data'])
    end
end

%start at time = 0
t0 = min(meas.t(1),inpt.t(1));
meas.t = meas.t - t0;
inpt.t = inpt.t - t0;

%compute steer angles


inpt.steer = state(:,3); %encoder ticks?

%compute wheel velocities (rad/s)
speed = state(:,2); %m/s
if opts.smooth
    %robust smooth, slow!
    speed = smooth(speed,5,'rlowess');
    
    set(figure,'name','speed')
    hold on
    plot(inpt.t,state(:,2),'r:')
    plot(inpt.t,speed,'b')
    plot(meas.t,meas.omit*max(speed),'k')
    xlabel('time (s)')
    ylabel('m/s')
    legend('raw','smooth','omit')
    makeLegible(14)
end

inpt.speed = speed;

inpt.rol=interp1(meas.t,meas.euler(:,1),inpt.t,'nearest','extrap');
inpt.pit=interp1(meas.t,meas.euler(:,2),inpt.t,'nearest','extrap');


%solve for initial yaw angle
if opts.inityaw
    ind0 = 1;
    indf = find(meas.s > 2,1);
    indf_yi = indf; %for plot
    
    resin={p,cinfo,mdl,meas,inpt,[ind0 indf]};

    [~,~,olog0]=residualRecbot(resin{:}); %for plot

    options = optimset('Display','iter');
    dyaw = fminbnd(@(dyaw) costyaw(resin,dyaw,sum(cinfo.inresidual)), -pi, pi, options);
    meas.euler(:,3)=meas.euler(:,3)+dyaw;
    
    resin{4}=meas;
    [~,~,olog]=residualRecbot(resin{:}); %for plot

end



if opts.plot    
    
    %%
    if 1
        %%
        %plot on aerial image
        
        zoom = 21;
        border = 1;
        [h_im,h_line] = plotPathOnAerial(meas.pos(:,1),meas.pos(:,2),meas.lat,meas.lon,zoom,border);
        plot(meas.pos(1,1),meas.pos(1,2),'b>')
        plot(meas.pos(end,1),meas.pos(end,2),'bs')
        
    end
    
    set(figure,'name','path')
    hold on
    h=[];
    h(end+1) = plot3(meas.pos(:,1), meas.pos(:,2), meas.pos(:,3),'DisplayName','meas');
    if opts.inityaw
        pos = olog0.pos_sensor;
        h(end+1) = plot3(pos(:,1),pos(:,2),pos(:,3),'r','DisplayName','pred, before yaw init');
        pos = olog.pos_sensor;
        h(end+1) = plot3(pos(:,1),pos(:,2),pos(:,3),'g','DisplayName','pred, after yaw init');
        plot3(pos(end,1),pos(end,2),pos(end,3),'gs')
        plot3(meas.pos(indf_yi,1),meas.pos(indf_yi,2),meas.pos(indf_yi,3),'bs')
    end
    xlabel('x')
    ylabel('y')
    zlabel('z')
    legend(h,'Location','Best')
    axis equal
    makeLegible(14)
    
    set(figure,'name','orientation')
    plot(meas.t,meas.euler*180/pi)
    legend('rol','pit','yaw')
    xlabel('time (s)')
    ylabel('deg')
    makeLegible(14)
    
    %to check synchronization, overlap
    set(figure,'name','speed')
    hold on
    meas_speed = numDiff(meas.s,meas.t);
    plot(meas.t,meas_speed,'b')
    plot(inpt.t,abs(inpt.speed),'r:')
    xlabel('time (s)')
    ylabel('m/s')
    legend('meas','inpt')
    makeLegible(14)
    
end

% blah = 1;

end

function c=costyaw(resin,dyaw,nr)

    %add yaw offset
    meas = resin{4};
    meas.euler(:,3) = meas.euler(:,3)+dyaw;
    resin{4} = meas;

    r = residualRecbot(resin{:});
    r = r(1:nr); %only systematic part
    c = r'*r;
end













