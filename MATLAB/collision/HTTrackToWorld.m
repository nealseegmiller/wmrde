function HT_track_to_world = HTTrackToWorld(mdl,HT_world)

HT_track_to_world = zeros(4,4,mdl.nt);
for tno = 1:mdl.nt
    fi = mdl.sprocketframeinds(tno); %frame index
    pfi = mdl.frames(fi).parent_ind; %parent frame index
    HT_track_to_world(:,:,tno) = HT_world(:,:,pfi)*mdl.frames(fi).HT_parent_jd0;
end