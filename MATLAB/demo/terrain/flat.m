function surface = flat()

%construct a flat surface, z=0 everywhere
surface = PlaneSurf([0 0 1 0]');
% surface=[]; %equivalent! update contact frames function interprets empty as z=0 plane
