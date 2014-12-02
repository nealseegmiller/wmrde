function [u,cmd] = rockyController(mdl, time, state)

%INPUTS
%OUTPUTS
%see template!

%control like an Ackerman steered car


%DIMENSIONS
k3 = 20/100;        %horizontal distance between D and wheels
k4 = 28.8/100;      %distance from D to steering axis of front wheels
k6 = 16/100;        %length of link from rocker joint to bogie joint
k9 = 139;           %angle of link between rocker and bogie joints (degrees)
k6x = sind(k9-90)*k6;

L = k4+k6x; %length (front to rear)
w = k3; %half track width

na = mdl.na;
nv = mdl.nf+5;

dlen = length(state)-mdl.nf;
%steering state indices (get steer angles)
% steer_si = namesToInds(mdl,{'S1','S2'}) + dlen;
steer_si = [3 9] + dlen;

radii = [mdl.frames(mdl.wheelframeinds).rad];


speed = 0.5;
turnrad = Inf;

% if time < 2
%     speed = 0.5;
%     turnrad = Inf;
% elseif time >= 2 && time < 6
%     speed = 0.5;
%     turnrad = -3;
% else
%     speed = 0.5;
%     turnrad = 3;
% end


if isinf(turnrad)
    omega = 0;
    gamma=zeros(2,1);
    vbl=speed;
    vbr=speed;
    vfl=speed;
    vfr=speed;
else
    omega = speed/turnrad;

    %steer angles [left right]
    gamma = atan([L/(turnrad-w) L/(turnrad+w)])';

    %steering limits
    lim = 60*pi/180;
    gamma(gamma > lim) = lim;
    gamma(gamma < -lim) = -lim;

    %wheel speeds (m/s)
    %TODO, CHECK THIS!
    vbl = speed - omega*w; %back left
    vbr = speed + omega*w; %back right
    vfl = (speed - omega*w)*cos(gamma(1)) + omega*L*sin(gamma(1)); %front left
    vfr = (speed + omega*w)*cos(gamma(1)) + omega*L*sin(gamma(1)); %front right
end

u = zeros(na,1);

%for indexing u
steer_ai = [1 5]; %steering rates
wheel_ai = [2:4,6:8]; %wheel vel
u(steer_ai) = 1/.2*(gamma - state(steer_si));


u(wheel_ai) = [vfl vbl vbl vfr vbr vbr]'./radii'; %wheel speeds (rad/s)

if nargout > 1
    cmd=zeros(nv,1); %commanded qvel
    cmd(4)=speed; %vx
    cmd(3)=omega; %wz
    cmd(mdl.actframeinds + 5) = u;
end






