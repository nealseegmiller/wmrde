% test_algebra.m
%test homogeneous transform and spatial vector algebra functions

%%
%test_linalg3
clear_all
close all
clc


A=[1 2 3; 4 5 6; 7 8 9]';
B=A+1;
disp('C=A*B')
C=A*B

%%
%test_transform
clear_all
close all
clc

%compose HT
euler=[10 20 30]'*pi/180;
translation=[1 2 3]';
HT = poseToHT(euler,translation)

%rotation matrices
disp('Rot = Rotz(yaw)*Roty(pit)*Rotx(rol)')
Rot = Rotz(euler(3))*Roty(euler(2))*Rotx(euler(1))

disp('inv(HT)=')
disp(inv(HT))

disp('HT*HT=')
disp(HT*HT)

disp('inv(HT)*HT=')
disp(inv(HT)*HT)


%apply HT
p=[2 3 4]'
disp('q=HT*p')
q=applyHT(HT,p)
disp('q=inv(HT)*p')
q=applyInvHT(HT,p)


vel=[1 2 3]'
[eulerrate,T_vel_to_eulerrate] = velToEulerrate(euler,vel)
disp('convert back:')
vel = eulerrateToVel(euler,eulerrate)


%test at singularity
disp('test at singularity:')
euler=[0 -pi/2 0]'
vel=[0 0 1]'
eulerrate = velToEulerrate(euler,vel)



%%
%test_spatial
clear_all
close all
clc


euler=[10 20 30]'*pi/180;
translation=[1 2 3]';
HT=poseToHT(euler,translation)

%HT to Plucker transform
P=HTToPlucker(HT)


%inv(HT) to Plucker
n=1e8;
%Matlab is slower than cpp, reduce the number of iterations so timing is feasible, & scale up elapsed time accordingly
irf=1e-4; %iteration reduction factor

if 0
    tic
    for i=1:(n*irf)
        P2 = invHTToPlucker(HT);
    end
    t=toc;
else
    %faster
    HT_ = repmat(HT,[1,1,n*irf]);
    tic
    P2_ = invHTToPlucker(HT_);
    t=toc;
    P2 = P2_(:,:,1);
end

disp('P2=inv(P)')
disp(P2)
fprintf('iterations: %e\n',n)
t_=t/irf*1e3;
fprintf('Elapsed time (ms): %f\n',t_)



m=[1 2 3 4 5 6]'
f=[2 3 4 5 6 7]'

%cross products

disp('crossMotion(m)*f=')
disp(crossMotion(m)*f)

disp('crossForce(m)*f')
disp(crossForce(m)*f)

%multiply spatial vector by Plucker transform
disp('P*m=')
disp(P*m)

disp('P''*m=')
disp(P'*m)


mass=10;
L=1;
W=.8;
H=.25;
cm=[L/2 0 H/2]';
I = inertiaBox(mass,L,W,H);
Is = toSpatialInertia(mass,cm,I)

n=1e7;
tic
for i=1:(n*irf)
    Is2=P'*Is*P;
end
t=toc;
disp('P''*Is*P=')
disp(Is2)
fprintf('iterations: %e\n',n)
t_=t/irf*1e3;
fprintf('Elapsed time (ms): %f\n',t_)


