function surface = ditch()


%for Crusher
lx = 30;
ly = 10;
ppm = 5;
size0 = 9;
ni = 4;
hmult = 4;
hmin = -Inf;
sig = [];
% sig = .1*ppm;

cx = 0; %center
cy = 0;

Z = fractalmat(lx,ly,ppm,size0,ni,hmult,hmin,sig);

%add the ditch
w = 8; %width
bw = 2; %basin width
d = 3; %depth

[X,~] = ndgrid(...
    linspace(cx-lx/2,cx+lx/2,size(Z,1)),...
    linspace(cy-ly/2,cy+ly/2,size(Z,2)));

%basin
Zd = zeros(size(X));
Zd(X >= -bw/2 & X <= bw/2) = -d;

%sides
m = 2*d/(w-bw);
b = -(d + m*bw/2);
I = X >= -w/2 & X < -bw/2;
Zd(I) = -m*X(I) + b;

I = X > bw/2 & X <= w/2;
Zd(I) = m*X(I) + b;

%make GridSurf
surface = GridSurf(...
    cx-lx/2,cx+lx/2,size(Z,1),...
    cy-ly/2,cy+ly/2,size(Z,2),Z+Zd);

end


