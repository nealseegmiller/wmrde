function [f,x,y] = rsgeng2D(N,rL,h,cl,opt)
%
% [f,x,y] = rsgeng2D(N,rL,h,cl,opt) 
%
% generates a square isotropic 2-dimensional Gaussian random rough surface 
% f(x,y) with NxN surface points. The surface has a Gaussian height distribution 
% function and Gaussian autocovariance functions (in both x and y), 
% where rL is the length of the surface side, h is the RMS height and cl is 
% the correlation length.
%
% Input:    N   - number of surface points (along square side)
%           rL  - length of surface (along square side)
%           h   - rms height
%           cl  - correlation length
%           opt - optional parameter (type 'plot' for plotting surface
%                 profile)
%
% Output:   f  - surface heights
%           x  - surface points
%           y  - surface points
%
% Last updated: 2009-02-23 (David Bergström).  
%

% format long;

x = linspace(-rL/2,rL/2,N); y = linspace(-rL/2,rL/2,N);

Z = h.*randn(N,N); % uncorrelated Gaussian random rough surface distribution
                   % with mean 0 and standard deviation h
                                   
% Gaussian filter
F = zeros(N,N); 
for n = 1:N;
  for m = 1:N;
      F(n,m) = 2/sqrt(pi)*exp(-((x(1,n))^2+(y(1,m))^2)/(cl^2/2));
  end;
end;

%correlation of surface including convolution (faltung), inverse
%Fourier transform and suitable prefactors

f = rL/N/cl*ifft2(fft2(Z).*fft2(F));

% optional plotting
if nargin<5 || isempty(opt)
    return;
end;
if nargin==5 
    if ischar(opt)
        if strcmp(opt,'plot');
            mesh(y,x,f);
            xlabel('x')
            ylabel('y')
            zlabel('z=f(x,y)')
            title('Plot of generated surface z=f(x,y)')
        else fprintf('%s is not a valid option. Type \''help rsgeng2D\'' for further details.\n',opt); 
        end;
    else fprintf('Option must be a string. Type \''help rsgeng2D\'' for further details.\n');
    end;
end;


