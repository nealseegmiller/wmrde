function v = triuv(M)
%copy the upper triangular portion of matrix to a vector
%by row!
%if M = [1 2 3; 4 5 6; 7 8 9], triuv(M) = [1 2 3 5 6 9]';
%compare to MATLAB triu()

M = M';
v = M(tril(true(size(M))));

