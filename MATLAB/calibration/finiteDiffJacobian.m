function [J,r,R,olog] = finiteDiffJacobian(fh,p,free,dp,addl_inputs)
%fh, function handle, vector output

if length(dp)==1
    dp = dp*ones(size(p));
end
[r,R,olog] = feval(fh,p,addl_inputs{:});

J = zeros(length(r),length(p));

for i=find(free);

    p_ = p;
    p_(i) = p_(i) + dp(i);
    r_ = feval(fh,p_, addl_inputs{:});
    J(:,i) = (r-r_)/dp(i);
end

%tried parfor, slower!