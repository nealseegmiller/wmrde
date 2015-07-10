%test wheel-ground contact models
clear_all
close all

ns = 51;
na = 51;
s = linspace(-.99,.99,ns);
alpha = linspace(-89,89,na)*pi/180;

[S,A] = ndgrid(s,alpha);
s = S(:)';
alpha = A(:)';

n=length(s);
dz = -.03*ones(1,n);
vz = zeros(1,n);



%%
wgc_fh = @ishigamiWgc; %function handle
wgc_name = func2str(wgc_fh);


%solve for Vx,Vy,Rw such that norm([Vx,Vy,Rw]) = 1;
%depends on definition of slip ratio
switch func2str(wgc_fh)
    case {'odeWgc','pacejkaWgc'}
        %s=(Rw-Vx)/Vx
        Vx = sqrt( 1 ./ ( (s+1).^2 + 1 + tan(alpha).^2 ) );
    case {'ishigamiWgc','ishigamiLUTWgc'}
        %s=(Rw-Vx)/Rw
        Vx = sqrt( (1-s).^2 ./ (1 + (1-s).^2.*(1+tan(alpha).^2) ) );
    otherwise
        disp('function handle not recognized')
end

Vy = Vx.*tan(alpha);
Rw = sqrt( 1 - Vx.^2 - Vy.^2 );

% Vx = -Vx;
% Rw = -Rw;

%contact point velocities:
vx = Vx-Rw;
vy = Vy;

Kp = 1e3;

p=wgcParams(wgc_fh,Kp);
feval(wgc_fh,p); %DEBUGGING, make lookup table

% return

[f,J]=feval(wgc_fh,p,[vx;vy;vz],Rw,dz);

%compute Jacobian numerically
if 0
    Jnum = zeros(size(J));
    myeps = 1e-6;
    for i = 1:5
        if i==1
            fp = feval(wgc_fh,p,[vx+myeps;vy;vz],Rw,dz);
            fm = feval(wgc_fh,p,[vx-myeps;vy;vz],Rw,dz);
        elseif i==2
            fp = feval(wgc_fh,p,[vx;vy+myeps;vz],Rw,dz);
            fm = feval(wgc_fh,p,[vx;vy-myeps;vz],Rw,dz);
        elseif i==3
            fp = feval(wgc_fh,p,[vx;vy;vz+myeps],Rw,dz);
            fm = feval(wgc_fh,p,[vx;vy;vz-myeps],Rw,dz);
        elseif i==4
            fp = feval(wgc_fh,p,[vx;vy;vz],Rw+myeps,dz);
            fm = feval(wgc_fh,p,[vx;vy;vz],Rw-myeps,dz);
        elseif i==5
            fp = feval(wgc_fh,p,[vx;vy;vz],Rw,dz+myeps);
            fm = feval(wgc_fh,p,[vx;vy;vz],Rw,dz-myeps);
        end      
        
%         Jnum(:,i,:) = reshape( (fp-f)/myeps, [3 1 n]); %forward diff
        Jnum(:,i,:) = reshape( (fp-fm)/(2*myeps), [3 1 n]); %central diff
    end
    
%     J = Jnum; %DEBUGGING
    J = J-Jnum; %DEBUGGING
end


        

Fx = reshape(f(1,:),ns,na);
Fy = reshape(f(2,:),ns,na);
Fz = reshape(f(3,:),ns,na);

%%
%PLOT


if 0
    set(figure, 'name', 'Vx')
    mesh(S,A*180/pi,reshape(Vx,ns,na))
    xlabel('s')
    ylabel('\alpha (deg)')
    title('V_x')

    set(figure, 'name', 'Vy')
    mesh(S,A*180/pi,reshape(Vy,ns,na))
    xlabel('s')
    ylabel('\alpha (deg)')
    title('V_y')

    set(figure, 'name', 'Rw')
    mesh(S,A*180/pi,reshape(Rw,ns,na))
    xlabel('s')
    ylabel('\alpha (deg)')
    title('R \omega')
end

%%
normalize = 1;
myaxis_ = [-1 1 -90 90]*1.05;

az = -45;
el = 30;

set(figure, 'name', 'Fx')
if normalize
    mesh(S,A*180/pi,Fx./Fz)
    zlabel('norm. lon. force')
else
    mesh(S,A*180/pi,Fx)
    zlabel('lon. force')
end
xlabel('slip ratio')
ylabel('slip angle (deg)')

axis tight
myaxis = axis;
myaxis(1:4) = myaxis_;
axis(myaxis);
makeLegible(14)
view(az,el)

if 0
    %%
    FileName = [wgc_name '_fx_vs_slip'];
    hgsave(gcf,[FileName '.fig'])
    export_fig([FileName '.pdf'])
end

set(figure, 'name', 'Fy')
if normalize
    mesh(S,A*180/pi,Fy./Fz)
    zlabel('norm. lat. force')
else
    mesh(S,A*180/pi,Fy)
    zlabel('lat. force')
end
xlabel('slip ratio')
ylabel('slip angle (deg)')

axis tight
myaxis = axis;
myaxis(1:4) = myaxis_;
axis(myaxis);
makeLegible(14)
view(az-90,el)

if 0
    %%
    FileName = [wgc_name '_fy_vs_slip'];
    hgsave(gcf,[FileName '.fig'])
    export_fig([FileName '.pdf'])
end


set(figure, 'name', 'Fz')
mesh(S,A*180/pi,Fz)
zlabel('normal force')

xlabel('slip ratio')
ylabel('slip angle (deg)')


axis tight
myaxis = axis;
myaxis(1:4) = myaxis_;
axis(myaxis);
makeLegible(14)
view(az+90,el)

if 0
    %%
    FileName = [wgc_name '_fz_vs_slip'];
    hgsave(gcf,[FileName '.fig'])
    export_fig([FileName '.pdf'])
end


%%
if 0
    %PLOT JACOBIANS
    dnamesrow={'dfx','dfy','dfz'};
    dnamescol={'dvx','dvy','dvz','dRw','ddz'};

    for i = 1:3 %row
        for j = 1:5 %col
            set(figure, 'name', [dnamesrow{i} '/' dnamescol{j}])
            mesh(S,A*180/pi,reshape(J(i,j,:),ns,na))
            xlabel('s')
            ylabel('\alpha (deg)')
%             pause(1)
        end
    end
end

