function C = covop(c)
%Covariance of the outer product of a vector of random variables
%e.g. for 3x1 random vector:
%[a b c]'*[a b c] = [aa ab ac; ba bb bc; ca cb cc];
%input c=cov([a b c]')
%output C=cov([aa ab ac bb bc cc]')

%Cov(ab,cd) = Cov(a,c)*Cov(b,d) + Cov(a,d)*Cov(b,c)
%***assumes (a,b,c,d) normally distributed and mean 0***

persistent m0 n I J

m = size(c,1);

if isempty(m0) || m~=m0
    n = m*(m+1)/2; %http://en.wikipedia.org/wiki/Triangular_number
    col=ones(m,1)*(1:m);
    row=col';
	I=triuv(row);
	J=triuv(col);
end

C = zeros(n);

for i = 1:n
    for j = i:n
        C(i,j) = c(I(i),I(j))*c(J(i),J(j)) + c(I(i),J(j))*c(J(i),I(j));
		C(j,i) = C(i,j); %symmetry
    end
end

m0=m; %for persistent