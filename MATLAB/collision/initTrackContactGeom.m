function contacts = initTrackContactGeom(frames)

%frames has fields .rad, .rad2, .L


%OPTIONS
npflat = 3;
% sides = [1 1 1 1] == 1;
% sides = [1 1 0 1] == 1;
sides = [1 0 0 0] == 1;


nt = length(frames); %number of tracks

contacts(1,nt) = TrackContactGeom();
for tno = 1:nt
    HT_contact_to_track = HTContactToTrack(frames(tno).rad, frames(tno).rad2, frames(tno).L, npflat, sides);
    contacts(tno).HT_track = HT_contact_to_track;
    
    contacts(tno).dvc_dvt = DvcontactDvtrack(HT_contact_to_track, frames(tno).rad);

end


if 0
    %%
    %DEBUGGING, PLOT
    tno = 1;
    HT_contact_to_track = contacts(tno).HT_track;
    
    np = size(HT_contact_to_track,3);
    cp = zeros(3,np);
    cp(:) = HT_contact_to_track(1:3,4,:);
    
    set(figure,'name',mfilename())
    hold on
    
    plot3(cp(1,:),cp(2,:),cp(3,:),'b.-','DisplayName','points')
    
    X = zeros(3,np);
    X(:) = HT_contact_to_track(1:3,1,:);
    h = drawQuivers(cp,X*.05);
    set(h,'Color','r','DisplayName','x axes')
    
    N = zeros(3,np);
    N(:) = HT_contact_to_track(1:3,3,:);
    h = drawQuivers(cp,N*.05);
    set(h,'Color','g','DisplayName','z axes')
    
    legend('Location','Best')
    
    axis equal
    axis tight
    
    view(30,30)
    d = .05;
    axis_lim = axis + d*[-1 1 -1 1 -1 1];
    axis(axis_lim)
    
    view(0,0)
    
    xlabel('x')
    ylabel('y')
    zlabel('z')

    makeLegible(14)
    

end

