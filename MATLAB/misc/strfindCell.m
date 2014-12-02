function I = strfindCell(text,pattern)
%INPUT
%text:      cell array of strings to search in
%pattern:   string, or cell array of strings to find
%OUTPUT
%I:         1 x size(pattern), text{I(i)} contains the substring pattern{i}
%           only returns the first occurence but warns if more than one

if ~iscell(pattern)
    pattern={pattern};
end

I=zeros(size(pattern));
for pno = 1:length(pattern)
    for tno = 1:length(text)
        k = strfind(text{tno},pattern{pno});
        if ~isempty(k) 
            if I(pno)==0
                I(pno) = tno;
            else
                disp([mfilename ': more than one occurence of pattern ' pattern{pno}])
                break
            end
        end
    end
end