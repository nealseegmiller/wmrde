function [f,J,fh] = uniformWgc(p,~,vc,Rw,dz)

%this function assumes the same wheel-ground contact model & parameters for all wheels,

f=[];
J=[];

fh=@odeWgc;
% fh=@pacejkaWgc;
% fh=@ishigamiWgc;
% fh=@ishigamiLUTWgc;

if nargin==0
    %just return function handle
    return
end

switch nargout
    case 0
        feval(fh,p); %initialize
    case 1
        f = feval(fh,p,vc,Rw,dz);
    case 2
        [f,J] = feval(fh,p,vc,Rw,dz);
end


