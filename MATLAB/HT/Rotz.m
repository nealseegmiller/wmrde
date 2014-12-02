function R = Rotz(angle)
    c = cos(angle);
    s = sin(angle);

    R = [c -s 0;
         s  c 0;
         0  0 1];
end