function Vq = changeRes(V,new)
%change the resolution of V by interpolation
%V:         m x n values array
%new:       1x2 size of output
%TODO, implement for other numbers of dimensions

F = griddedInterpolant(V,'linear','none');
old = size(V);

if length(new) == 2
    x1q = linspace(1,old(1),new(1));
    x2q = linspace(1,old(2),new(2));
    [X1q,X2q] = ndgrid(x1q,x2q);
    Vq = F(X1q(:),X2q(:));
    Vq = reshape(Vq(:),new(1),new(2));
else
    disp([mfilename ': not implemented for dim' num2str(length(new))])
end

