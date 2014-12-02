function surface = plane()

%creates an oriented plane surface

euler = [-5 0 0]*pi/180; %Euler angles
R = euler_to_rot(euler);

pt = [0 0 0]';

surface = PlaneSurf(R(:,3),pt);



