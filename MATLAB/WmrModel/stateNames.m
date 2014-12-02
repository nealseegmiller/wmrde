function out = stateNames(mdl)
%mdl:   WmrModel object
%out:   1xns cell array of strings, names of each element

olen = SIZEORIENT();
dlen = -1 + olen + 3;
ns = mdl.nf + dlen;

[isorient,ispos,isjd] = stateIs(ns);
out = cell(1,ns);
if olen==3
    out(isorient) = {'rol','pit','yaw'};
elseif olen==4
    out(isorient) = {'q1','q2','q3','q4'};
end
out(ispos) = {'x','y','z'};
out(isjd) = {mdl.frames(2:end).name};
