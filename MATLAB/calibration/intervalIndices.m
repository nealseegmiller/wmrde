function [ind0,indf] = intervalIndices(s,last_ind0,d_ind0,d_indf,skip,omit)
%INPUT
%s:         n x 1, distance or time
%last_ind0: ind0 on the previous iteration
%d_ind0:    distance/time between segments (>0 required)
%d_indf:    segment length/duration (>0 required)
%skip:      n x 1 logical, do not start or end on these indices
%omit:      n x 1 logical, interval must not contain these indices

n=length(s);

ind0=n;
for i = last_ind0+1:n-1
    if (last_ind0 == 0 || s(i)-s(last_ind0) >= d_ind0) && ~skip(i) && ~omit(i)
        ind0=i;
        break
    end
end

indf=n;
for i = ind0+1:n-1
    if (s(i)-s(ind0) >= d_indf && ~skip(i)) || omit(i+1)
        indf=i;
        break
    end
end