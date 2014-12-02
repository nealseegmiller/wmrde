function R = orientToRot(orient)

if numel(orient)==3
    R = eulerToRot(orient);
elseif numel(orient)==4
    R = quatToRot(orient);
end