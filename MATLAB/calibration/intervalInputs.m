function [time,inpt] = intervalInputs(t1,tf,dt,inpt)


fields = fieldnames(inpt);
I = find(~strcmp(fields,'t'))';

%interpolate inputs
time = (t1:dt:tf)';
time(end) = tf;

%loop over all fields of inpt (except .t) and interpolate
for i = I
    inpt.(fields{i}) = interp1(inpt.t,inpt.(fields{i}),time,'nearest','extrap');
end
inpt.t = time;

