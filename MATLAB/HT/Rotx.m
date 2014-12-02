function R = Rotx(angle)
    c = cos(angle);
    s = sin(angle);

    R = [1 0  0; 
         0 c -s;
         0 s  c];
end