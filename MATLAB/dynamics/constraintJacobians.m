function [A, cis] = constraintJacobians(mdl, state0, qvel0, vis_act, HT_world, contacts)
%Compute constraint Jacobians for forward Dynamics

%PRE-PROCESS

%get from WmrModel object
ns = length(state0);
nv = mdl.nf+5; %number of elements in qvel
[~,~,isjd] = stateIs(ns);

%contacts
incontact = [contacts.incontact];
npic = sum(incontact); %number of points in contact

%actuators
na = sum(vis_act); %number of actuated DOF

%END PRE-PROCESS

%contact constraints
ncc = 0;
if npic > 0
    if mdl.nw > 0
        Acontact = wheelJacobians(mdl, HT_world, contacts);
    elseif mdl.nt > 0
        Acontact = trackJacobians(mdl, HT_world, contacts);
    end    
    ncc=3*npic; %number of contact constraints
end

%actuator constraints
nac = na; %number of actuator constraints
Aact = eye(nv); 
Aact = Aact(vis_act,:);

%joint constraints
njc=0;
if ~isempty(mdl.hjc_fh)
    [~,Jc] = feval(mdl.hjc_fh, mdl, state0(isjd), qvel0(7:end));
    
    njc = size(Jc,1); %number of holonomic joint constraints
    Ajoint = [zeros(njc,6), Jc]; %Jacobian wrt qvel
end

nc = ncc + nac + njc; %total number of constraints

%allocate
A = zeros(nc,nv); %A*qacc=b
cis = struct; %logicals for indexing constraints

%wheel constraints
cis.contact = false(1,nc); %is wheel constraint
if ncc > 0
    cis.contact(1:ncc) = true; 
    A(cis.contact,:) = Acontact;
end

%actuator constraints
cis.act = false(1,nc); %is actuator constraint
if nac > 0
    cis.act((ncc+1):(ncc+nac)) = true; 
    A(cis.act,:) = Aact;
end

%holonomic joint constraints
cis.joint = false(1,nc); %is joint constraint
if njc > 0;
    cis.joint((ncc+nac+1):end)=true;
    A(cis.joint,:) = Ajoint;
end

