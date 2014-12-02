function out = isValidHT(HT)
    %check if HT is a valid homogenous transform
    out = true;

    if any(size(HT) ~= 4)
        disp([mfilename ': must be a 4x4 matrix'])
        out = false;
    end
    
    %check that rotation matrix is valid
    if ~isValidRot(HT(1:3,1:3))
        disp([mfilename ': rotation matrix invalid'])
        out = false;
    end

    if any(HT(4,1:4) ~= [0 0 0 1]);
        disp([mfilename ': last row must be [0 0 0 1]'])
        out = false;
    end

end