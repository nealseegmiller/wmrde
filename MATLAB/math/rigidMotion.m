function [R,t]=rigidMotion(a,b,origin_a,origin_b)
%Solve for optimal rotation (R) and translation (t) between 2D or 3D points using SVD
% http://nghiaho.com/?page_id=671
%INPUT
%a,b:                   2xn or 3xn matrices of points
%origin_a,origin_b:     (optional) 2x1 or 3x1 vectors. If provided solve only for rotation.
%OUTPUT
%R,t:                   Rotation, translation. b = R*a + t;

n=size(a,2); %number of points

if nargin < 3
    %make origins the centroids
    origin_a = mean(a,1);
    origin_b = mean(b,1);
end

ac=a-origin_a*ones(1,n);
bc=b-origin_b*ones(1,n);

H=ac*bc';
[U,~,V] = svd(H);
R=V*U';

if det(R) < 0
    R(:,end)=-R(:,end);
end

t=origin_b-origin_a;