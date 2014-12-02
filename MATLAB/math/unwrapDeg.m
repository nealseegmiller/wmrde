function out = unwrapDeg(in)
%INPUTS
%in, vector of wrapped angles

out = zeros(size(in));
out(1) = in(1);
for i = 2:numel(in)
    out(i) = out(i-1) + diffDeg(in(i),in(i-1));
end
    
    