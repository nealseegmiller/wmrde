function opts = simopts()

opts.dt = .05; %time step size
opts.nsteps = 400; %number of time steps
opts.model_fh = @zoe; %or other function in /demo/model
opts.ideal_actuators = true; %if true, actuators can supply any force required to achieve desired speed
opts.terrain_fh = @terrain;  %you can modify the terrain function in /demo/terrain
opts.init_in_contact = true; %if true, init suspension joint positions such that wheels are in contact with terrain.
opts.log = true; %must log simulation data for plots
opts.animate = true;
opts.use_builtin_solver = false;
opts.dyn = true; %if false, do kinematic simulation only
opts.profile = false;
opts.plot = true;