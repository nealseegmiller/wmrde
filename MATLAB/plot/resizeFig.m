function resizeFig(xs,ys,h_fig)
%resize figure window without changing position
%xs:        x scale (1 for unchanged, 2 to double size, etc.)
%ys:        y scale
%h_fig:     (optional) figure handle. if not provided resize current figure

if nargin < 3
    h_fig = gcf;
end

mypos = get(h_fig,'Position');
dx = mypos(3)*(xs-1);
dy = mypos(4)*(ys-1);
mypos([1 3]) = mypos([1 3]) + dx*[-1 1];
mypos([2 4]) = mypos([2 4]) + dy*[-1 1];
set(h_fig,'Position',mypos);