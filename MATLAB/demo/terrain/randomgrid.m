function surface = randomgrid()

% cx = 0; 
% cy = 0;

% %for Crusher
% lx = 30;
% ly = 10;
% ppm = 5;
% RmsHeight = .2;
% CorrLength = 1;

% %for Zoe
lx = 20;
ly = 10;
ppm = 10;
RmsHeight = .15;
CorrLength = .8;

cx = 8;
cy = 0;

% for LandTamer
% lx = 40;
% ly = 20;
% ppm = 10;
% RmsHeight = .10;
% CorrLength = .8;

% cx = 16; %center
% cy = 7;

% %for Rocky7
% lx = 10;
% ly = 5;
% ppm = 20;
% RmsHeight = .1;
% CorrLength = 0.8;

Z = randommat(lx,ly,ppm,RmsHeight,CorrLength);

surface = GridSurf(...
    cx-lx/2,cx+lx/2,size(Z,1),...
    cy-ly/2,cy+ly/2,size(Z,2),Z);


