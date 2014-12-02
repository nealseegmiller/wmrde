function Z = randommat(lx,ly,ppm,RmsHeight,CorrLength)
%INPUTS
%lx:            length in x (m)
%ly:            length in y (m)
%ppm:           points per meter
%RmsHeight:     RMS height (m)
%CorrLength:    correlation length (m)

SurfLength = max(lx,ly);
N = ppm*SurfLength+1;

Z = rsgeng2D(N,SurfLength,RmsHeight,CorrLength,[]);

%crop if necessary
if lx ~= ly
    Nx = lx*ppm+1;
    Ny = ly*ppm+1;
    Z = Z(1:Nx,1:Ny);
end
