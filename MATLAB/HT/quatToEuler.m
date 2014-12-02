function euler = quatToEuler(q)

euler = rotToEuler(quatToRot(q));

%TODO, a direct way?