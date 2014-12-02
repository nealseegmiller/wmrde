function S = skew3(v)
%turns a 3x1 vector into a skew symmetric matrix (for computing cross
%products)
%a x b = skew3(a)*b;

S = [    0, -v(3),  v(2)
      v(3),     0, -v(1)
     -v(2),  v(1),    0];