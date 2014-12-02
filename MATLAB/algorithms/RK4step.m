function [y,out]=RK4step(odefh,h,time,y0,addlin)

[k1,out]=feval(odefh,time,y0,addlin{:});
[k2,~]=feval(odefh,time+h/2,y0+h/2*k1,addlin{:});
[k3,~]=feval(odefh,time+h/2,y0+h/2*k2,addlin{:});
[k4,~]=feval(odefh,time+h,y0+h*k3,addlin{:});

y=y0+h/6*(k1 + 2*k2 + 2*k3 + k4);
