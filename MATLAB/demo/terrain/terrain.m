function surfaces = terrain()
%output a cell array of surface objects

surfaces = {};

% surfaces{end+1} = flat();
% surfaces{end+1} = plane();
% surfaces{end+1} = ramp();
surfaces{end+1} = randomgrid();
% surfaces{end+1} = fractalgrid();
% surfaces{end+1} = ditch();
% surfaces = [surfaces, rocks(surfaces)];