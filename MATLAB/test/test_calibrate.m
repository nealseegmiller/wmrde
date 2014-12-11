% test_calibrate.m

clear_all
close all
clc

%OPTIONS
opts.loadmatdata = 1; %to save time, load meas, inpt from .mat file instead of reading .txt & preprocessing
opts.save_cinfo = 0; %save CalibrationInfo
opts.inittosaved = 0;

opts.dyn = 0; %dynamic sim (else kinematic)

opts.d_indf = 10; %segment length
opts.d_ind0 = 10; %distance between segments

opts.dp = 1e-4; %delta p for Jacobian, 
opts.md_tol = Inf;

opts.uncal = 1; %compute residuals using uncalibrated parameter values

opts.printresults = 1;
opts.plot = 1;
opts.plot_scatter = 1;
opts.plot_yaw_hist = 1;
opts.plot_params = 1;

%%

if 1
    %LandTamer
    opts.model_fh = @landtamer;
    opts.initcal_fh = @initcalSkidsteer;
    opts.readdata_fh = @readdataLandtamer;
    opts.residual_fh = @residualSkidsteer;
    
%     testname = 'Taylor_may13_30_dirt'; opts.starttime = 1053961201;
%     testname = 'Taylor_may13_30_grass'; opts.starttime = 1053964566;
    testname = 'Taylor_may13_30_parkinglot'; opts.starttime = 1053966982;

    opts.stoptime = opts.starttime + 180;
    
    ldir = [datalogdir() 'LandTamer\' testname '\'];
    opts.filenames = {
        [ldir testname '_commands.data'], ...
        [ldir testname '_pose.data'], ...
        [ldir testname '_state.data']
    };
else
    %RecBot
    opts.model_fh = @recbot;
    opts.initcal_fh = @initcalRecbot;
    opts.readdata_fh = @readdataRecbot;
    opts.residual_fh = @residualRecbot;
    
    testname = 'Taylor_may13_asphalt'; opts.starttime = 416030;
%     testname = 'Taylor_may13_dirt'; opts.starttime = 417182; %missing first ~1200 seconds of state data?
%     testname = 'Taylor_may13_grass'; opts.starttime = 421495; %missing first ~550 seconds of state data?

    opts.stoptime = opts.starttime + 180;
    
    ldir = [datalogdir() 'RecBot\' testname '\'];
    opts.filenames = {
        [ldir testname '_novatel.data'], ...
        [ldir testname '_status.data']
    };
end

%make WmrModel object
mdl = feval(opts.model_fh);
ns = mdl.nf-1+6; %number of elements in state, assume Euler angles


%initialize calibration
[cinfo, x_unc, P_unc, Q] = feval(opts.initcal_fh,mdl,opts.dyn);

np = length(x_unc); %number of parameters
nr = sum(cinfo.inresidual); %size of residual (not including stochastic)

x = x_unc;
P = P_unc;

if opts.inittosaved
    saved = load('_autosave/calibinfo.mat');
    x(:) = saved.x;
    P(:) = saved.P;
end

%load data
%meas must have n x 1 fields: .t, .s, .skip, .omit
if opts.loadmatdata
    load('_autosave/logdata.mat') %meas, inpt, logname
else
    %x or x0?
    [meas, inpt] = feval(opts.readdata_fh,opts.filenames,opts.starttime,opts.stoptime,x_unc,cinfo,mdl);
    save('_autosave/logdata.mat','meas','inpt')
end
n = length(meas.t);

%indicate bad indices
if ~isfield(meas,'skip')
    meas.skip = false(n,1); 
end
if ~isfield(meas,'omit')
    meas.omit = false(n,1);
end

%DATA LOGGING
nres = nr + (nr^2+nr)/2; %size of residual (including stochastic)
clog = CalibrationLog(np,nres,1);

set_md_tol(clog,opts.md_tol);
set_x0(clog,x);
set_P0(clog,P);

if opts.uncal
    clog_uncal = CalibrationLog(np,nres,1);
end


%initialize
s = meas.s;
% s = meas.t; %split into intervals by time vs. distance
d_ind0_ = 0;

%additional inputs for residual function
ind0=0;
indf=0;
addl_inputs = {cinfo, mdl, meas, inpt, [ind0 indf], opts.dyn};
I_idx = 5; %index of I in addl_inputs

iter = 0;

free = ~cinfo.isfixed; 

while indf < n

    [ind0,indf] = intervalIndices(s,ind0,d_ind0_,opts.d_indf,meas.skip,meas.omit);
    
    if ind0 == indf
        continue
    end
    
    iter = iter+1;
        
    d_ind0_=opts.d_ind0;
    
    I = [ind0 indf];
    addl_inputs{I_idx} = I;
    
    %process update (add noise)
    ds = meas.s(indf) - meas.s(ind0);
    P = P + Q*ds^2;
    
    %measurement update
    [J,res,C,olog] = finiteDiffJacobian(opts.residual_fh,x,free,opts.dp,addl_inputs);

    if ~any(cinfo.isstoch(free))
        %no free stochastic parameters, so ignore stochastic residual
        J(nr+1:end,:) = 0;
    end
    [x,P,md] = kfMeasUpdate(x,P,res,J,C,opts.md_tol);
    
    
    if any(x(cinfo.isstoch) < 0)
        disp([mfilename ': making stochastic parameters positive']) %DEBUGGING
        x(cinfo.isstoch) = abs(x(cinfo.isstoch)); 
    end
    
    
    %DATA LOGGING
    set_I(clog,iter,I);
    set_res(clog,iter,res);
    set_C(clog,iter,C);
    set_odomlog(clog,iter,olog);
    set_x(clog,iter,x);
    set_P(clog,iter,P);
    set_mdist(clog,iter,md); %DEBUGGING, Mahalanobis distance
    
    
    if opts.uncal
        %compute residual using uncalibrated parameter values
        
        [res,C,olog] = feval(opts.residual_fh,x_unc,addl_inputs{:});
        set_res(clog_uncal,iter,res);
        set_C(clog_uncal,iter,C);
        set_odomlog(clog_uncal,iter,olog);
    end
    
    fprintf('processed inds %d to %d of %d\n',ind0,indf,n) %DEBUGGING
    
end

if opts.save_cinfo
    %%
    save('_autosave/calibinfo.mat','x','P','cinfo')
end

%POST-PROCESSING

%%

inliers = [];
% inliers = clog.mdist < opts.md_tol; %DEBUGGING

inityaw0 = true;
clog;

res_pos_cal = res_pos(clog,cinfo.inresidual,inityaw0,inliers);
res_yaw_cal = res_yaw(clog,cinfo.inresidual,inliers);

C_pos_cal = C_pos(clog,cinfo.inresidual,inityaw0,inliers);
C_yaw_cal = C_yaw(clog,cinfo.inresidual,inliers);

if opts.uncal
    res_pos_unc = res_pos(clog_uncal,cinfo.inresidual,inityaw0,inliers);
    res_yaw_unc = res_yaw(clog_uncal,cinfo.inresidual,inliers);
    
    C_pos_unc = C_pos(clog_uncal,cinfo.inresidual,inityaw0,inliers);
    C_yaw_unc = C_yaw(clog_uncal,cinfo.inresidual,inliers);
end

if opts.printresults && opts.uncal
    %%
    %Euclidean distance position error
    dist_cal = sqrt(sum(res_pos_cal.^2,2)); %calibrated
    dist_unc = sqrt(sum(res_pos_unc.^2,2)); %uncalibrated
    
    fprintf('segment length: %.1f, distance btwn segments: %.1f\n',opts.d_indf,opts.d_ind0)
    M = [mean(res_pos_unc(:,1)) var(res_pos_unc(:,1)) mean(res_pos_unc(:,2)) var(res_pos_unc(:,2)) mean(res_yaw_unc) var(res_yaw_unc) mean(dist_unc) mean(abs(res_yaw_unc));
         mean(res_pos_cal(:,1)) var(res_pos_cal(:,1)) mean(res_pos_cal(:,2)) var(res_pos_cal(:,2)) mean(res_yaw_cal) var(res_yaw_cal) mean(dist_cal) mean(abs(res_yaw_cal))];
    M(3,:) = (abs(M(1,:)) - abs(M(2,:)))./abs(M(1,:));
    disptable(M, {'mean x','var x','mean y','var y','mean yaw','var yaw','mean pos','mean |yaw|'},{'uncalibrated','calibrated','% improvement'}, '%.4f',2)
    
end




if opts.plot

    qNames = stateNames(mdl);

    
    if opts.plot_scatter
        %%
        %scatter plots of position residuals
        
        %axis limits
        tmp = res_pos_cal;
        if opts.uncal
            tmp = [tmp; res_pos_unc];
        end
        lim = max(abs(tmp)); 
        lim(1:2) = max(lim(1:2)); %same x,y limits
        lim = lim + .01;
        lim = [-lim(1) lim(1) -lim(2) lim(2) -lim(3) lim(3)]*1.1;
        
        set(figure,'name','scatter (calibrated)')
        hold on
        xlabel('x')
        ylabel('y')
        zlabel('z')
        
        conf = .95;
        dim = 3;
        [h,h_cent,he_line,he_surf] = plotScatter(res_pos_cal,conf,dim,gca);
        
        clear h_leg
        h_leg(1) = he_line(1);
        set(h_leg(1),'DisplayName',sprintf('observed cov, %.0f%%',conf*100))
        
        set(he_line,'Color','r')
        set(he_surf,'FaceColor','r')
        
        [he_line,he_surf] = drawEllipse(mean(C_pos_cal,3),zeros(3,1),36*2,conf,[]);
        
        h_leg(2) = he_line(1);
        set(h_leg(2),'DisplayName',sprintf('pred cov, %.0f%%',conf*100))
        
        set(he_line,'Color','g','LineStyle','--')
        set(he_surf,'FaceColor','g')

        legend(h_leg,'Location','Best')
        axis equal
        axis(lim);
        makeLegible(14)
        
        if opts.uncal
            set(figure,'name','scatter (uncalibrated)')
            hold on
            xlabel('x')
            ylabel('y')
            zlabel('z')
            
            [h,h_cent,he_line,he_surf] = plotScatter(res_pos_unc,conf,dim,gca);
            h_leg = he_line(1);
            set(h_leg,'DisplayName',sprintf('observed cov, %.0f%%',conf*100))
            
            set(he_line,'Color','r')
            set(he_surf,'FaceColor','r')
            
            legend(h_leg,'Location','Best')
            axis equal
            axis(lim);
            makeLegible(14)
            
        end
        


    end
    
    if opts.plot_yaw_hist
        %%
        
        set(figure,'name','yaw error hist')
        hold on
        ylabel('count')
        
        %bins
        tmp = res_yaw_cal;
        if opts.uncal
            tmp = [tmp; res_yaw_unc];
        end
        nbins = 30;
%         bins = nbins;
        bins = linspace(min(tmp*180/pi),max(tmp*180/pi),nbins);
        
        conf = .95;
        h_leg = [];

        if opts.uncal
            [~,~,h_hist] = plotHist(res_yaw_unc*180/pi,bins,[],gca);
            h_leg(end+1) = h_hist; set(h_leg(end),'DisplayName','uncalibrated')
            
        end
        [N,~,h_hist,h_bars] = plotHist(res_yaw_cal*180/pi,bins,conf,gca);
        set(h_hist,'FaceColor','g')
        set(h_bars,'Color','g','LineStyle','--')
        
        h_leg(end+1) = h_hist; set(h_leg(end),'DisplayName','calibrated')
        h_leg(end+1) = h_bars(1); set(h_leg(end),'DisplayName',sprintf('observed %.0f%% err',conf*100))
        
        
        %predicted yaw variance
        clear h_bars
        sig = sqrt(chi2inv(conf,1));
        stdev = mean(sqrt(C_yaw_cal)*180/pi);
        h_bars(1) = plot(-sig*stdev*[1 1], [0 max(N)]);
        h_bars(2) = plot( sig*stdev*[1 1], [0 max(N)]);
        
        set(h_bars,'Color','g','LineWidth',1.5)
        
        h_leg(end+1) = h_bars(1); set(h_leg(end),'DisplayName',sprintf('pred %.0f%% err',conf*100))
        
        legend(h_leg,'Location','Best')
        
        makeLegible(14)
        

    end

    
    if opts.plot_params && any(~cinfo.isfixed)
        %%
        set(figure, 'Name', 'param values')
        hold on
        plot(0:1:iter, clog.x0x(:,~cinfo.isfixed),'.-')
        legend(cinfo.param_names(~cinfo.isfixed),'Location','Best')
        xlabel('iteration')
        makeLegible(14)
        
        set(figure, 'Name', 'param variance')
        hold on
        plot(0:1:iter, clog.P0P_diag(:,~cinfo.isfixed),'.-')
        legend(cinfo.param_names(~cinfo.isfixed),'Location','Best')
        xlabel('iteration')
        makeLegible(14)
        
    end
    
    
end
