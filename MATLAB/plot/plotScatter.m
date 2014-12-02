function [h,h_cent,he_line,he_surf] = plotScatter(res,conf,dim,h_axis)
%plot 3D residual scatter
%INPUTS
%r:         nx3 matrix of residuals [x,y,z]
%(optional)
%conf:      confidence for ellipse
%dim:       {2,3} for 2D or 3D plot
%h_axis:    existing axis handle
%OUTPUTS


if nargin < 4
    h_axis = [];
    if nargin < 3
        dim = 2;        
        if nargin < 2
            conf = [];
        end
    end
end

if isempty(h_axis)
    set(figure,'name','scatter')
    xlabel('x')
    ylabel('y')
    zlabel('z')
else
    axes(h_axis)
end

m = mean(res,1); %centroid

if dim == 2
    h = plot(res(:,1),res(:,2),'.');
    h_cent = plot(m(1),m(2),'+');
    marg = 3;
elseif dim == 3
    h = plot3(res(:,1),res(:,2),res(:,3),'.');
    h_cent = plot3(m(1),m(2),m(3),'+');
    marg = [];
end

if ~isempty(conf)
    %ellipse for scatter covariance
    nvertices=36*2;

    C = res'*res/size(res,1); %scatter covariance, assumes mean = 0
    [he_line,he_surf]=drawEllipse(C,zeros(3,1),nvertices,conf,marg);
end

if isempty(h_axis)
    axis equal
    make_legible(14)
end




