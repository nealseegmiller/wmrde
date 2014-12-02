%generate optimized c code to transform spatial inertia, P'*I*P;

clear_all
close all
clc

%matrix blocks
% M = [M1 M3 
%      M2 M4]


%spatial inertia matrix
%block 1 is symmetric
I1 = sym('I1_%d%d',[3 3]); I1 = sym(I1,'real');
I1(2,1) = I1(1,2);
I1(3,1) = I1(1,3);
I1(3,2) = I1(2,3);


I2 = sym('I2_%d%d',[3 3]); I2 = sym(I2,'real');
%block 2 is skew-symmetric
I2(1,1) = 0;
I2(2,2) = 0;
I2(3,3) = 0;
% why slower?
% I2(2,1) = -I2(1,2);
% I2(3,1) = -I2(1,3);
% I2(3,2) = -I2(2,3);

%block 3 is the transpose of block 2
I3 = I2';

I4 = sym('I4_%d%d',[3 3]); I4 = sym(I4,'real');
%block 4 is diagonal
tmp = I4(1,1);
I4(:) = 0;
I4(1,1) = tmp;
I4(2,2) = tmp;
I4(3,3) = tmp;
I = [I1,I3; I2,I4];


%Plucker transform
P1 = sym('P1_%d%d',[3 3]); P1 = sym(P1,'real');
P2 = sym('P2_%d%d',[3 3]); P2 = sym(P2,'real');
P3 = sym(zeros(3));
P4 = P1;
P = [P1,P3; P2,P4];

ccode((P'*I*P), 'file', 'temp.txt');

% return

%%
clc

% inds = [
%  0, 3, 6, 18,21,24
%  1, 4, 7, 19,22,25
%  2, 5, 8, 20,23,26
% 
%  9,12,15, 27,30,33
% 10,13,16, 28,31,34
% 11,14,17, 29,32,35
% ];

%with extra elements for alignment
inds = [
 0, 4, 8, 24,28,32
 1, 5, 9, 25,29,33
 2, 6,10, 26,30,34

12,16,20, 36,40,44
13,17,21, 37,41,45
14,18,22, 38,42,46
];

%get from inputs
sym_names = {};
for i=1:numel(I)
    sym_name = char(I(i));
    if I(i) ~= 0 && ~any(strcmp(sym_name,sym_names)) && sym_name(1) ~= '-'
        disp(['Real ' sym_name ' = I[' num2str(inds(i)) '];'])
        sym_names{end+1} = sym_name; %#ok<SAGROW>
    end
end

fprintf('\n')

sym_names = {};
for i=1:numel(P)
    sym_name = char(P(i));
    if P(i) ~= 0 && ~any(strcmp(sym_name,sym_names))
        disp(['Real ' sym_name ' = P[' num2str(inds(i)) '];'])
        sym_names{end+1} = sym_name; %#ok<SAGROW>
    end
end

fprintf('\n')

temp = fileread('temp.txt');
temp = strrep(temp,char(10),'');

% declare temporary variables
% how many?
ti = strfind(temp,'t');
tno = zeros(size(ti));
for i=1:length(ti)
    tno(i) = sscanf(temp(ti(i)+1:end),'%d');
end
nt = max(tno); %number of temp vars

fprintf('Real ')
for i = 2:nt
    if rem(i,10) == 0
        fprintf('\n'); 
    end
    fprintf('t%d,',i)
end
fprintf('\b;\n\n')

%fix output
k=0;
for i = 0:1:5
    for j = 0:1:5
        temp = strrep(temp,sprintf('A0[%d][%d]',i,j),sprintf('R[%d]',inds(i+1,j+1)));
    end
end

disp(temp)
    
    





