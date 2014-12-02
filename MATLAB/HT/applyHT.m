function q = applyHT(HT,p)

    %q = HT*p
    %p must be a 3xn matrix!

    R = HT(1:3,1:3);
    t = HT(1:3,4);
    
    t = t*ones(1,size(p,2));
    q = R*p + t;

end