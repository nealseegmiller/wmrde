function surfs = terrain()
%output a cell array of surface objects

surfs = {};

% surfs{end+1} = flat();
% surfs{end+1} = plane();
% surfs{end+1} = ramp();

surfs{end+1} = randomgrid();
printGridSurf(surfs{end},[resourcedir() 'gridsurfdata.txt'])


% surfs{end+1} = fractalgrid();
% surfs{end+1} = ditch();
% surfs = [surfs, rocks(surfs)];