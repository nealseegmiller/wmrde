function printGridSurf(s,filename)
%print GridSurf data to txt file
%s:         GridSurf object
%filename:  string

Zdata = s.Z.Values;
fmt = '%.9f,';
fid = fopen(filename,'w');

fprintf(fid,[fmt fmt '%d,'],s.X(1,1),s.X(end,1),size(Zdata,1));
fprintf(fid,[fmt fmt '%d,'],s.Y(1,1),s.Y(1,end),size(Zdata,2));
for i = 1:numel(Zdata)
   fprintf(fid,fmt,Zdata(i));
end

fclose(fid);

return