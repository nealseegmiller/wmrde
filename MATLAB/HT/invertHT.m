function ht = invertHT(ht)
    %inv(ht) may be faster?
    
    R = ht(1:3,1:3);
    t = ht(1:3,4);
    
    ht(1:3,1:3) = R';
    ht(1:3,4) = -R'*t;
    
end