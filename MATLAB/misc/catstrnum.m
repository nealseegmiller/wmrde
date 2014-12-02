function out=catstrnum(str,num)
%out is a (1 x num) cell array of strings
%cell i is str concatenated with i
out=cell(1,num);
for i=1:num
    out{i}=[str,num2str(i)];
end
