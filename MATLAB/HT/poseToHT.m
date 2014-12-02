function HT = poseToHT(orient,pos)
%get HT matrix from pose (orienation & position)

if length(orient) == 3
    R = eulerToRot(orient);
elseif length(orient) == 4
    R = quatToRot(orient);
end

HT = eye(4);
HT(1:3,1:3) = R;
HT(1:3,4) = pos;

