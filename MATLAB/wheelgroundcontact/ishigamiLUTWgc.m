function [f,J] = ishigamiLUTWgc(p,vc,Rw,dz)
%LUT- lookup table f = f(s,beta,dz)
%only 3 calls to interp3

persistent FX FY FZ

%range limits
s_lim = 1;
beta_lim = 89.9*pi/180; %singularity at beta=0
dz_llim = -.08;
dz_ulim = 0;

if nargin==1
    
    if 1
        %ignore p - load from memory
        load('_autosave/LUT.mat')
        return
    end
    
    %make the lookup tables given the parameters
    ns = 101;
    nb = 101;
    nz = 51;
    n = ns*nb*nz;
    
    s = linspace(-s_lim,s_lim,ns); 
    beta = linspace(-beta_lim,beta_lim,nb); 
    dz = linspace(dz_llim,dz_ulim,nz);
    
    [S,B,DZ] = ndgrid(s,beta,dz);
    s = S(:)';
    beta = B(:)';
    dz = DZ(:)';
    
    %solve for Vx,Vy,Rw such that norm([Vx,Vy,Rw]) = 1
    %forward not reverse
    %s = (Rw-Vx)/Rw
    Vx = sqrt( (1-s).^2 ./ (1 + (1-s).^2.*(1+tan(beta).^2) ) );
    Vy = Vx.*tan(beta);
    Rw = sqrt( 1 - Vx.^2 - Vy.^2 );
    
    vx = Vx-Rw;
    vy = Vy;
    vz = zeros(1,n);
    
    f = ishigamiWgc(p,[vx;vy;vz],Rw,dz);
    
    FX = reshape(f(1,:),ns,nb,nz);
    FY = reshape(f(2,:),ns,nb,nz);
    FZ = reshape(f(3,:),ns,nb,nz);
    
    FX = griddedInterpolant(S,B,DZ,FX,'linear');
    FY = griddedInterpolant(S,B,DZ,FY,'linear');
    FZ = griddedInterpolant(S,B,DZ,FZ,'linear');
    
    VarNames = {'S','B','DZ','FX','FY','FZ'};
    save('_autosave/LUT.mat','p',VarNames{:})
    
    if 0
        %%
        %PLOT
        M={FX.Values, FY.Values, FZ.Values};

        FigNames = {'FX','FY','FZ'};
        for i=1:3
            j= round(.8*nz);
            set(figure,'name',FigNames{i})
            mesh(S(:,:,j),B(:,:,j)*180/pi,M{i}(:,:,j)) 
            xlabel('s')
            ylabel('beta (deg)')
            makeLegible(14)
        end
    end
    
    if 0
        %%
        %print to txt file
        FileNames = {'LUT_FX_zoe','LUT_FY_zoe','LUT_FZ_zoe'};
%         FileNames = {'LUT_FX_rocky','LUT_FY_rocky','LUT_FZ_rocky'};
        
        for fileno = 1:3
            F = f(fileno,:);
            
            fmt = '%.9f,';
            fid = fopen([resourcedir() FileNames{fileno} '.txt'],'w');

            fprintf(fid,[fmt fmt '%d,'],s(1),s(end),ns);
            fprintf(fid,[fmt fmt '%d,'],beta(1),beta(end),nb);
            fprintf(fid,[fmt fmt '%d,'],dz(1),dz(end),nz);
            for i = 1:numel(F)
               fprintf(fid,fmt,F(i));
            end

            fclose(fid);
        end
        
    end
    blah = 1;
    return
end



vx=vc(1,:);
vy=vc(2,:);
vz=vc(3,:);

if nargout > 1
    %append for Jacobian
    %finite difference for derivatives wrt vx,vy,Rw,dz
    nw = length(vx);
    
    myeps = 1e-4;
    
    vx = [vx, vx+myeps, vx, vx, vx];
    vy = [vy, vy, vy+myeps, vy, vy];
    vz = [vz, vz, vz, vz, vz]; %compute analytically
    Rw = [Rw, Rw, Rw, Rw+myeps, Rw];
    dz = [dz, dz, dz, dz, dz+myeps];
    
end

   
Kd = p(1); %damping

%convert to slip ratio and angle
[s,beta] = calcSlip(vx,vy,Rw,2);

%stay in bounds to avoid NaN during interpolation
s(s>s_lim) = s_lim;
s(s<-s_lim) = -s_lim;
beta(beta>beta_lim) = beta_lim;
beta(beta<-beta_lim) = -beta_lim;
dz(dz>dz_ulim) = dz_ulim;
dz(dz<dz_llim) = dz_llim;

%use griddedInterpolant class
fx=FX(s,beta,dz);
fy=FY(s,beta,dz);
fz=FZ(s,beta,dz);

%COPIED FROM NON-LUT FUNCTION
%add damping to normal force
fz = fz - Kd*vz;

%normal force must be positive
fz(fz<0)=0;

%reverse sign for reverse, TODO check this!
Vx=vx+Rw;
isrev = Rw<0;
fx(isrev) = -fx(isrev); 
isrev = Vx<0;
fy(isrev) = -fy(isrev);

f=[fx;fy;fz];

if nargout > 1
    f=f(:,1:nw);
    
    J = zeros(3,5,nw);
    
    i=0;
    for j = [1 2 4 5] %loop over col
        i=i+1;
        %forward difference
        I = i*nw + (1:nw);
        
        J(1,j,:) = (fx(I)-fx(1:nw))/myeps;
        J(2,j,:) = (fy(I)-fy(1:nw))/myeps;
        J(3,j,:) = (fz(I)-fz(1:nw))/myeps;
    end
  
    J(3,3,:) = -Kd; %wrt vz
    J(3,:,fz(1:nw)==0)=0;

end

return


