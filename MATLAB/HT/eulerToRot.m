function R = eulerToRot(euler)
%euler = [rol pit yaw]'
%ASSUME 321 ROTATION ORDER! THIS IS THE SAE CONVENTION

rol = euler(1);
pit = euler(2);
yaw = euler(3);

if 0
    R = Rotz(yaw)*Roty(pit)*Rotx(rol);
else
    %faster. derived using symbolic toolbox
    sR=sin(rol); cR=cos(rol);
    sP=sin(pit); cP=cos(pit);
    sY=sin(yaw); cY=cos(yaw);
    R=[...
        cP*cY, cY*sP*sR - cR*sY, sR*sY + cR*cY*sP
        cP*sY, cR*cY + sP*sR*sY, cR*sP*sY - cY*sR
          -sP,            cP*sR,            cP*cR
    ];
end