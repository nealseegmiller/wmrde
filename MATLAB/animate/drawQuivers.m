function h = drawQuivers(xyz,uvw,h)
%xyz:   3xn matrix, position
%uvw:   3xn matrix, vector specifying direction

if nargin > 2
    %handle exists, just update 'XData','YData,'ZData'
    set(h, 'XData', xyz(1,:), 'YData', xyz(2,:), 'ZData', xyz(3,:), ...
        'UData', uvw(1,:), 'VData', uvw(2,:), 'WData', uvw(3,:) )
else
    %handle does not exist, plot from scratch
    h = quiver3(xyz(1,:), xyz(2,:), xyz(3,:), uvw(1,:), uvw(2,:), uvw(3,:), ...
        'AutoScale', 'off', 'MaxHeadSize', .5);
end