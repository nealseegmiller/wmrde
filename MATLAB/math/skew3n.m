function S = skew3n(V)

%turns n 3x1 vectors into n skew symmetric matrices. see skew3
%V is 3xn
%S is 3x3xn

n = size(V,2);

S = zeros(9,n);
S(2,:) = V(3,:);
S(3,:) = -V(2,:);
S(4,:) = -V(3,:);
S(6,:) = V(1,:);
S(7,:) = V(2,:);
S(8,:) = -V(1,:);

S=reshape(S,[3 3 n]);
