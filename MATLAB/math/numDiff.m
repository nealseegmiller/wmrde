function xdot = numDiff(x,t,method)
%numerically differentiate with respect to time (or something else)
%angle values should be unwrapped prior to using this function
%INPUTS
%x, nxm
%t, nx1 or [], time
%method, 
%   1 forward difference (default)
%   0 central difference
%  -1 backward difference
%OUTPUTS
%xdot, nxm d/dt x

if nargin < 3
    method = 1;
    if nargin < 2
        t = [];
    end
end

if method==1
    %forward difference
    csx=circshift(x,-1);
    dx = csx - x;

    if isempty(t)
        xdot = dx;
    else
        dt = circshift(t,-1) - t;
        xdot = dx./( dt*ones(1,size(x,2)) );
    end

    xdot(end,:) = xdot(end-1,:);
    
elseif method==0
    %central difference
    csx_m=circshift(x,-1); %minus
    csx_p=circshift(x,1); %plus

    dx = csx_m - csx_p;
    
    if isempty(t)
        xdot = dx;
    else
        dt = circshift(t,-1) - circshift(t,1);
        xdot = dx./( dt*ones(1,size(x,2)) );
    end
    
    xdot(1,:) = xdot(2,:);
    xdot(end,:) = xdot(end-1,:);
    
elseif method==-1
    %backward difference
    csx=circshift(x,1);
    dx = x - csx;
    
    if isempty(t)
        xdot = dx;
    else
        dt = t - circshift(t,1);
        xdot = dx./( dt*ones(1,size(x,2)) );
    end

    xdot(1,:) = xdot(2,:);
    
else
    disp([mfilename ': method invalid'])
    
end

