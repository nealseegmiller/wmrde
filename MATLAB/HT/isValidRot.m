function out = isValidRot(R)
    %check if rotation matrix is valid
    
    out = true;
    
    if any(size(R) ~= 3)
        out = false;
    end
    
    %determinant == 1 and orthogonal to within tolerance
    tol = 1e-10;
    err = abs(R'*R - eye(3));
    if abs(det(R) - 1) > tol || any(err(:) > tol)
        out = false;
    end
end