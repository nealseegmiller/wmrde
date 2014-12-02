function contacts = updateModelContactGeom(mdl, surfaces, HT_world, min_npic, contacts)
%update all ContactGeom objects for a WMR model
%subfunction used by both odeKin and odeDyn
%INPUTS
%min_npic:  minimum number of points in contact. mdl.min_npic for kinematic sim, 0 for dynamic


if mdl.nw > 0
    radii = [mdl.frames(mdl.wheelframeinds).rad];    
    contacts = updateWheelContactGeom(surfaces,HT_world(:,:,mdl.wheelframeinds),radii);
elseif mdl.nt > 0
    HT_track_to_world = HTTrackToWorld(mdl,HT_world);
    contacts = updateTrackContactGeom(surfaces,HT_track_to_world,contacts);
end
contacts = whichPointsInContact(contacts, min_npic);