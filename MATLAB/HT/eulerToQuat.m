function q = eulerToQuat(euler)

q = rotToQuat(eulerToRot(euler));

%TODO, a direct way?