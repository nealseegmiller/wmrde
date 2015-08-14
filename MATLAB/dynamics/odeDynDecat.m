function [s,q,e]= odeDynDecat(y,nf,na)
%decatenate input vector y (or output vector ydot) of the odeDyn function
%INPUTS
%y:     y (or ydot)
%nf:    number of frames in WmrModel object
%na:    number of actuated joints
%OUTPUTS
%s:     state (or statedot)
%q:     qvel (or qacc)
%e:     interr (or err)

nv=nf+5; %length(q)
ns=size(y,1)-nv-na; %length(s)

s = y(1:ns,:);
q = y(ns+1:ns+nv,:);
e = y(ns+nv+1:end,:);
