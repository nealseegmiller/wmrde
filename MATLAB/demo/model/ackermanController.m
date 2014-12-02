function [u,cmd] = ackermanController(mdl, time, state)
%controller for rear wheel drive Ackerman steered vehicle (Gator)
%INPUTS
%OUTPUTS
%see template!

%control like an Ackerman steered car

persistent L wf wr na avi steer_si
if isempty(L)
    %DIMENSIONS
    wfi = mdl.wheelframeinds; %[fl fr bl br]
    
    steer_fi = [mdl.frames(wfi(1:2)).parent_ind]; %steering frame indices
    tmp = steer_fi;
    susp_fi = [mdl.frames(steer_fi).parent_ind]; %suspension frame indices
    if all(susp_fi > 1)
        tmp = susp_fi;
    end
        
    L =  mdl.frames(tmp(1)).HT_parent_jd0(1,4); %wheelbase
    wf = mdl.frames(tmp(1)).HT_parent_jd0(2,4); %1/2 track width, front
    wr = mdl.frames(wfi(3)).HT_parent_jd0(2,4); %1/2 track width, rear
    
    avi = mdl.actframeinds + 5; %actuated velocity indices (in qvel)
    
    dlen = length(state)-mdl.nf;
    steer_si = steer_fi + dlen; %steering state indices (steer angles)
    
end

radii = [mdl.frames(wfi).rad];

if time < 1
    speed = .5;
    turnrad = Inf;
elseif time >= 1 && time < 5
    speed = 1;
    turnrad = Inf;
elseif time >= 5 && time < 10
    speed = 1;
    turnrad = -5;
else
    speed = 1;
    turnrad = 5;
end


if isinf(turnrad)
    omega = 0;
    gamma=zeros(2,1);
    vbl=speed;
    vbr=speed;
else
    omega = speed/turnrad;

    %steer angles [left right]
    gamma = atan([L/(turnrad-wf) L/(turnrad+wf)])';

    %steering limits
    lim = 60*pi/180;
    gamma(gamma > lim) = lim;
    gamma(gamma < -lim) = -lim;

    %wheel speeds (m/s)
    vbl = omega*(turnrad-wr); %back left
    vbr = omega*(turnrad+wr); %back right

end


na=length(avi); %number of actuated joints
u = zeros(na,1);
%assume steer frames before wheel frames
u(1:2) = 1/.2*(gamma - state(steer_si)); %steering rates
u(3:4) = [vbl vbr]' ./ radii(3:4)'; %wheel vel (rad/s)

if nargout > 1
    cmd=zeros(mdl.nf+5,1); %commanded qvel
    cmd(4)=speed; %vx
    cmd(3)=omega; %wz
    cmd(avi)=u;
end

