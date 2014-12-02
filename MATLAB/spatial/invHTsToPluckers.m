function out = invHTsToPluckers(in)
%converts the inverse of n 4x4 homogeneous transforms to 6x6 Plucker transforms

%in 4x4xn matrix
%out 6x6xn matrix

%in(:,:,i) = [E, -E*r; 0,  1]

n = size(in,3);

out = zeros(6,6,n);

Et = permute(in(1:3,1:3,:),[2 1 3]); %E' permute transposes 3x3 matrices
mEr = zeros(3,n); mEr(:) = in(1:3,4,:); %-E*r

out(1:3,1:3,:) = Et;
out(4:6,4:6,:) = Et;
temp = skew3n(mEr);
for i = 1:n
    temp(:,:,i) = -Et(:,:,i)*temp(:,:,i); %TODO, possible to vectorize this?
end
out(4:6,1:3,:) = temp;
