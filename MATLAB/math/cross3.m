function c = cross3(a,b)
%computes the cross product of two length 3 vectors
%returns a 3x1 vector
%c = a x b

% if numel(a) ~= 3 || numel(b)  ~= 3
%     disp('(cross3.m) inputs must be length 3 vectors')
%     return
% end

c = zeros(3,1);
c(1) = -a(3)*b(2) + a(2)*b(3);
c(2) = a(3)*b(1) - a(1)*b(3);
c(3) = -a(2)*b(1) + a(1)*b(2);

