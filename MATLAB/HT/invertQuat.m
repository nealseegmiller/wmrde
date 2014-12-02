function invq = invertQuat(q)
%compute inverse (or reciprocal) of quaternion
%http://en.wikipedia.org/wiki/Quaternion

invq = [q(1); -q(2:4)];
invq = invq/sum(q.^2); %normalize


