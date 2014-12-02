function out = PluckerToHT(in)

%converts a 6x6 Plucker transform to a 4x4 homogeneous transform

E = in(1:3,1:3);
mErx = in(4:6,1:3);			% - E r cross
out = [ E, unskew3(mErx*E'); 0 0 0 1 ];