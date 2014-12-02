function out = jointSpaceVelNames(mdl)
%mdl:   WmrModel object
%out:   1xns cell array of strings, names of each element

out = {'wx','wy','wz','vx','vy','vz',mdl.frames(2:end).name};