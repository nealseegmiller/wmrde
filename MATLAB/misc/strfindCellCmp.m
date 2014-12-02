function I = strfindCellCmp(text,pattern)
%INPUT
%text:      cell array of strings to search in
%pattern:   string, or cell array of strings to find
%OUTPUT
%I:         1 x size(pattern), text{I(i)} is identical to pattern{i} (using strcmp)
%           only returns the first occurence but warns if more than one

if ~iscell(pattern)
    pattern={pattern};
end

I = zeros(size(pattern));
for i=1:length(pattern)
    ind = find(strcmp(text,pattern{i})); %must be identical
    if ~isempty(ind)
        I(i) = ind(1);
    end
    if numel(ind) > 1
        disp([mfilename ': more than one occurence of pattern ' pattern{i}])
    end
end

