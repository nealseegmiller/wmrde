function surface = fractalgrid()

%for Crusher
% lx = 30;
% ly = 10;
% ppm = 5;
% size0 = 9;
% ni = 4;
% hmult = 5;
% hmin = -.25;


%for Zoe
lx = 20;
ly = 10;
ppm = 10;
size0 = 9;
ni = 4;
hmult = 3;
hmin = -.15;


cx = 4; %center
cy = 0;

Z = fractalmat(lx,ly,ppm,size0,ni,hmult,hmin);

surface = GridSurf(...
    cx-lx/2,cx+lx/2,size(Z,1),...
    cy-ly/2,cy+ly/2,size(Z,2),Z);


