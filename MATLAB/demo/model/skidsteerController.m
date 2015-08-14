function [u,cmd] = skidsteerController(mdl, time, ~)
%INPUTS
%OUTPUTS
%see template!


%TODO, make track width a WmrModel property?
%track width, reference y location of front left wheel
B = 2*abs( mdl.frames(2).HT_parent_jd0(2,4) ); %Crusher

if mdl.nw > 0
    %wheeled
    nw = mdl.nw;
    radii = [mdl.frames(mdl.wheelframeinds).rad];
elseif mdl.nt > 0
    %tracked
    nw = mdl.nt; 
    radii = [mdl.frames(mdl.sprocketframeinds).rad];
end

%no actuated joints other than wheels
nv = mdl.nf+5;


%for yaw instability
% speed = 5;
% omega = 2;
%   
% if time < 1
%     speed = .1*speed;
%     omega = .1*omega;
% elseif time >= 1 && time < 3
% 
% elseif time >= 3
%     speed = .1*speed;
%     omega = .1*omega;
% end

%for rollover
% if time < 1
%     speed = 5;
%     omega = 0;
% elseif time >= 1
%     speed = 5;
%     omega = 7;
% end

% for ramp test
if time < 1
    speed = 1;
    omega = 0;
elseif time >= 1 && time < 10
    speed = 2;
    omega = 0;
elseif time >= 10 && time < 15
    speed = 2;
    omega = 1;
else
    speed = 2;
    omega = -1;
end

%for monte-carlo sim.
% if time < 5
%     speed = 2;
%     omega = .5;
% else
%     speed = 2;
%     omega = -.5;
% end

% speed = .5;
% omega = 0;

% [vl vr]' = [.5 .5; -1/B 1/B]\[speed omega]';

vl = [1 -B/2]*[speed omega]';
vr = [1  B/2]*[speed omega]';

u = zeros(nw,1);
u(1:2:end) = vl; %odd numbered wheels are on left
u(2:2:end) = vr; %even numbered wheels are on right
u = u./radii';


cmd=zeros(nv,1); %commanded qvel
cmd(4)=speed;
cmd(3)=omega;
cmd(mdl.actframeinds + 5) = u;









