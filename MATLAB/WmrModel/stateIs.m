function [isorient,ispos,isjd] = stateIs(ns)
%logicals to index state vector
%INPUTS
%ns:        size of state vector
%OUTPUTS 
%all 1 x ns logicals that specify which elements of state vector are:
%isorient:  body frame orientation (Euler angles or quaternion) 
%ispos:     body frame position
%isjd:      joint displacement

olen = SIZEORIENT();
plen = olen + 3; %pose length

isorient = false(1,ns);
isorient(1:olen) = true;
% isorient(4:plen)=true; %DEBUGGING

if nargout > 1
    ispos = false(1,ns);
    ispos(olen+1:plen) = true;
%     ispos(1:3)=true; %DEBUGGING

    if nargout > 2
        isjd = false(1,ns);
        isjd(plen+1:end) = true;
    end
end
