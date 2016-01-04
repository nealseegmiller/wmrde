function surfs = terrain()
%output a cell array of surface objects

%make terrain
surfs = {};

% surfs{end+1} = flat();
% surfs{end+1} = ramp(); %must also uncomment flat

%GridSurf

setseed(123);
surfs{end+1} = randomgrid();
% surfs{end+1} = fractalgrid();
% surfs{end+1} = ditch();
% printGridSurf(surfs{end},[resourcedir() 'gridsurfdata.txt'])



% surfs = [surfs, rocks(surfs)];