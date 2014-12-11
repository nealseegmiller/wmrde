function [u,cmd] = zoeController(mdl, time, state)

%INPUTS
%OUTPUTS
%see template!


persistent nv steer_fi wheel_vi
if isempty(nv)
    nv = mdl.nf+5;
    steer_fi = namesToInds(mdl,{'s1','s2'}); %indices of steer frames in WmrModel
    wheel_vi = mdl.wheelframeinds+5; %indices of wheel velocities in qvel
end
dlen = length(state)-mdl.nf; %5 if using Euler angles, 6 if using quaternions


%TODO, wheel acceleration limits?

speed = .5;
turnrad = 1000;
    
% if time < .5
%     speed = .5;
%     turnrad = 1000;
% elseif time >= .5 && time < 1
%     speed = 1;
%     turnrad = 1000;
% elseif time >= 1 && time < 5
%     speed = 1;
%     turnrad = 5;
% elseif time >= 5
%     speed = 1;
%     turnrad = -5;
% end

steer = state(steer_fi+dlen); %steer angles
[~,yawrate,wvel,steer_rates] = zoeArcController(speed,turnrad,steer);

u = wvel;

cmd=zeros(nv,1); %commanded qvel
cmd(4)=speed;
cmd(3)=yawrate;
cmd(wheel_vi)=u;
cmd(steer_fi+5) = steer_rates;

end




