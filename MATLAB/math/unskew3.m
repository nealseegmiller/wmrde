function v = unskew3(S)
%turns a 3x3 skew symmetric matrix back into a 3x1 vector
%see skew3

v = [S(6) -S(3) S(2)]';

