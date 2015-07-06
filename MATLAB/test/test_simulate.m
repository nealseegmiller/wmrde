%test_simulate.m

clear_all
close all
clc

%OPTIONS
opts.dyn = 1; %dynamic sim (else kinematic)
%for dynamic sim
opts.use_erp_cfm = 0; %use error reduction parameters & constraint force mixing like Open Dynamics Engine
opts.ideal_actuators = 0;

opts.initincontact = 1;
opts.log = 1;
opts.animate = 1;
opts.plot = 0;
opts.profile = 0;

dt = .04;
nsteps = 10/dt + 1;

%uncomment one of the following:
% model_fh = @zoe;
model_fh = @rocky;
% model_fh = @talon;

%make WmrModel object
if opts.animate
    [mdl, state, qvel, anim] = feval(model_fh);
else
    [mdl, state, qvel] = feval(model_fh);
end

feval(mdl.wgc_fh,mdl.wgc_p); %DEBUGGING, initialize wheel-ground contact model

mdl.bsm_fh = []; %DEBUGGING
if opts.use_erp_cfm
    mdl.wgc_fh = [];
end
if opts.ideal_actuators
    mdl.act_fh = [];
end

nf = mdl.nf;
nw = mdl.nw;
na = mdl.na;

%make terrain
surfs = {};

% surfs{end+1} = flat();
% surfs{end+1} = ramp(); %must also uncomment flat

setseed(123);
surfs{end+1} = randomgrid();
% surfs{end+1} = fractalgrid();
printGridSurf(surfs{end},[resourcedir() 'gridsurfdata.txt'])



%init contact geometry
if mdl.nw > 0
    contacts(1:mdl.nw) = WheelContactGeom();
    radii = [mdl.frames(mdl.wheelframeinds).rad];
elseif mdl.nt > 0
    contacts = initTrackContactGeom(mdl.frames(mdl.sprocketframeinds));
end
np = sum([contacts.np]);

%initialize wheel/terrain contact
if opts.initincontact
    state = initTerrainContact(mdl,surfs,contacts,state);
end

state0 = state

interr = zeros(na,1);
time = 0;

if opts.log
    dlog = SimulationLog(nsteps,nf,np,na,opts.dyn);

    %log initial values
    set_time(dlog,1,time);
    set_state(dlog,1,state);
    set_qvel(dlog,1,qvel);
end


if opts.animate
    figure(anim.h_fig)
    resizeFig(1.5,1.5)
    
    drawSurfaces(surfs,anim.h_axis);
    
    makeLegible(14)
    axis off
    
%     view(0,0)
    view(30,30)
    
    zoom(1.5)
    
end

do_pause = true;

%simulate

if opts.profile
    profile on
end
for i = 2:nsteps
    
    if opts.dyn
        y=[state; qvel; interr];
        [ydot,out] = odeDyn(time,y,mdl,surfs,contacts,dt);
        y = y + ydot*dt;
        [state,qvel,interr]=odeDynSplitVec(y,nf,na);
    else
        [statedot,out] = odeKin(time,state,mdl,surfs,contacts);
        state = state + statedot*dt;
    end
    time = time + dt;
    

    if opts.log
        set_time(dlog,i,time);
        set_state(dlog,i,state);
        set_qvel(dlog,i,qvel);
        %see SimulationLog for more to log
    end

    if opts.animate
        updateAnimation(anim, out.HT_parent, out.contacts);
        
        if do_pause
            pause
            do_pause = false;
        end
    end

end
if opts.profile
    profile off
    profile viewer
end

state


if opts.plot
    %%
    olen = SIZEORIENT();
    [isorient, ispos] = stateIs(nf);
    qnames = stateNames(mdl);
    lw = 1.5; %line width
    
    
    %orientation vs time
    set(figure,'name','orientation vs time')
    hold on
    if olen==3
        plot(dlog.time, dlog.orient*180/pi,'LineWidth',lw)
        ylabel('deg')
    else
        plot(dlog.time, dlog.orient)
    end
    xlabel('time (s)')
    legend(qnames{isorient},'Location','Best')
    makeLegible(14)
    
end



