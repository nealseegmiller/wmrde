function contacts = whichPointsInContact(contacts,min_npic)
%sets .incontact field for WheelContactGeom or TrackContactGeom objects
%INPUTS
%contacts:  1xnw array of WheelContactGeom objects, or
%           1xnt array of TrackContactGeom objects
%min_npic:  scalar, minimum number of points in contact

%shortcut (for dynamic sim)
if min_npic == 0
    for wtno = 1:length(contacts)
        contacts(wtno).incontact = contacts(wtno).dz < 0;
    end
    return
end

tol = 1e3;

%get from objects
dz = [contacts.dz];
HT_contact_to_world = cat(3,contacts.HT_world);

if any(isnan(dz))
    disp([mfilename ': NaN error'])
    return
end

np = length(dz); %number of points
incontact = dz <= 0;
npic = sum(incontact); %number of points in contact

nwt = length(contacts); %number of wheels/tracks
%TODO, speed up?
whichwt = uncatIndexArray([contacts.np]); %which wheel/track each contact point belongs to


%shortcut. return if all points in contact
if npic == np
    for wtno = 1:nwt
        contacts(wtno).incontact = incontact(whichwt == wtno);
    end
    return
end 

% min_npic = sum(~isinf(dz)); %DEBUGGING! all points in contact

ind0 = max(npic,min_npic); %at least min_npic wheels must be in contact
[~,isort] = sort(dz);

for ind = ind0:np
    
    %update incontact
    incontact = false(1,np);
    incontact(isort(1:ind)) = true;
    npic = sum(incontact);
    
    %shortcut. break if all wheels in contact
    if npic == np
        break
    end

    %check contact points are stable
    cp = zeros(2,npic);
    cp(:) = HT_contact_to_world(1:2,4,incontact);
    if isStable(cp,tol)
        break
    end

end

for wtno = 1:nwt
    contacts(wtno).incontact = incontact(whichwt == wtno);
end

end

function out = isStable(cp,tol)
    %determine if stable via condition number of covariance of contact points
    npic = size(cp,2);
    
    cp_ = cp - mean(cp,2)*ones(1,npic); %subtract mean
    cp_cov = cp_*cp_'/(npic-1);
    ev = eigenvalues2x2(cp_cov); %compare to MATLAB eig()
    
    ev = abs(ev); %DEBUGGING, numerical issues
    val = ev(1)/ev(2);
%     disp(num2str(val)) %DEBUGGING
    
    if val < tol
        out = true;
    else
        out = false;
    end
end

function ev = eigenvalues2x2(M)
    %http://www.math.harvard.edu/archive/21b_fall_04/exhibits/2dmatrices/index.html
    T = M(1,1) + M(2,2); %trace
    D = M(1,1)*M(2,2) - M(1,2)*M(2,1); %determinant
    
    t1 = T/2;
    t2 = sqrt(T*T/4-D);
    
    ev(1) = t1 + t2;
    ev(2) = t1 - t2;
end



