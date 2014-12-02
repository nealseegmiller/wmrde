function out = structNoDependent(in) %#ok<STOUT>
%like MATLAB built in struct() function, converts object to struct
%but, to save memory, does not copy Dependent parameters

mc = metaclass(in);
for i = 1:length(mc.PropertyList)
    if ~mc.PropertyList(i).Dependent
        name = mc.PropertyList(i).Name;
        eval(['out.' name '= in.' name ';'])
    end
end