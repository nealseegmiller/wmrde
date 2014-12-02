% clear_all.m
%clear all except for breakpoints

%store breakpoints
tmp = dbstatus;
save('mybreakpoints.mat','tmp')

%clear all
clear all %clear all

%reload breakpoints
load('mybreakpoints.mat')
dbstop(tmp)

clear tmp
delete('mybreakpoints.mat')