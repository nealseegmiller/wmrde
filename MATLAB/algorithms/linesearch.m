function [x,cost,grad,alpha] = linesearch(x,cost,grad,p,alpha_max,fh_calc_cost,fh_calc_gradient)
%line search algorithm to satisfy strong Wolfe conditions
%http://en.wikipedia.org/wiki/Wolfe_conditions
%http://pages.cs.wisc.edu/~ferris/cs730/chap3.pdf
%http://terminus.sdsu.edu/SDSU/Math693a_s2004/Lectures/06/lecture-static-color.pdf
%http://sentientdesigns.net/math/mathbooks/Number%20theory/Numerical%20Optimization%20-%20J.%20Nocedal,%20S.%20Wright.pdf

%INPUTS
%x, nx1, cost function inputs
%cost, scalar cost
%grad, 1xn, gradient of cost wrt x
%p, nx1, search direction
%alpha_max, the maximum step size (recommend 1 for Newton's method)
%fh_calc_cost, fh_calc_gradient are function handles to nested functions, so they can access externally scoped variables
%calc_cost takes input x
%OUTPUTS
%new x = x + alpha*p
%cost, grad computed for new x
%alpha (between 0 and alpha_max)

%OPTIONS
alpha_tol = 1e-3; %required resolution for alpha in zoom
%for Wolfe conditions:
c1 = 1e-4;
c2 = 0.9;

%back up at alpha=0
x0 = x;
cost_0 = cost; %cost at alpha=0;
grada_0 = grad*p; %gradient wrt alpha at alpha=0, should always be negative!

%first try taking the maximum step
alpha_j = alpha_max;
x = x0 + alpha_j*p;
cost = feval(fh_calc_cost,x);
grad = feval(fh_calc_gradient);

cost_j = cost;
grada_j = grad*p;

dozoom = false;
if cost_j > cost_0 + c1*alpha_j*grada_0
    %decrease condition not satisfied
    dozoom = true;
    alpha_hi = alpha_j;
    alpha_lo = 0; 
    %grad, addlout already backed up
    cost_lo = cost_0;
else

    if abs(grada_j) <= c2*abs(grada_0) %Ferris
%     if abs(grada_j) <= -c2*grada_0 %Blomgren/sdsu
        %strong Wolfe conditions satisfied!
        alpha=alpha_j;
    else
        %curvature condition not satisfied
        if grada_j >= 0
            dozoom = true;
            alpha_hi = 0;
            alpha_lo = alpha_j;
            cost_lo = cost_j;
        else
            %strong Wolfe conditions not satisfied, step could be larger
            alpha=alpha_j;
        end
    end
end

if dozoom
%     disp('zoom') %DEBUGGING
    while true

        %PROBLEM, curvature condition may be impossible to satisfy if there is a discontinuity in the gradient at the minimum.
        %need a way to break
        if abs(alpha_hi-alpha_lo) < alpha_tol
            
            %close enough, break with alpha_lo
            if alpha_j ~= alpha_lo
                %need to call cost and gradient functions again, not just to set cost & grad
                %but also to set externally scoped vars required by parent function
                x = x0 + alpha_lo*p;
                cost = feval(fh_calc_cost,x);
                grad = feval(fh_calc_gradient);
            end
            alpha = alpha_lo;
            
            break
        end

        alpha_j = (alpha_lo + alpha_hi)/2; %bisection
        x = x0 + alpha_j*p;
        cost = feval(fh_calc_cost,x);
        cost_j = cost;

%         fprintf('alpha_j: %f, cost_j: %f, cost_lo: %f\n', alpha_j, cost_j, cost_lo)

        if cost_j > cost_0 + c1*alpha_j*grada_0 || cost_j >= cost_lo
            alpha_hi = alpha_j;
        else
            grad = feval(fh_calc_gradient);
            grada_j = grad*p;
            if abs(grada_j) <= -c2*grada_0 %Ferris, Blomgren
                %strong Wolfe conditions satisfied!
                alpha=alpha_j;
                break
            end
            if grada_j*(alpha_hi-alpha_lo) >= 0
                alpha_hi = alpha_lo;
            end
            alpha_lo = alpha_j;
            cost_lo = cost_j;
        end
    end
end



