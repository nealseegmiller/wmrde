function [qacc, u, out] = forwardDynUnc(mdl, state0, qvel0, u, HT_world, contacts, H, C)
%Compute *unconstrained* forward dynamics, add forces directly through tau

%TODO, ELIMINATE DUPLICATION
%PRE-PROCESS

%get from WmrModel object
ns = length(state0);
nv = mdl.nf+5; %number of elements in qvel
[~,~,isjd] = stateIs(ns);

%to compute wheel/sprocket radius * angular rate
if mdl.nw > 0
    wsframeinds = mdl.wheelframeinds; %wheel/sprocket frame indices
    whichwt = 1:mdl.nw; %which wheel/track does contact point belong to
elseif mdl.nt > 0
    wsframeinds = mdl.sprocketframeinds;
    whichwt = uncatIndexArray([contacts.np]);
end
radii = [mdl.frames(wsframeinds).rad];

%contacts
incontact = [contacts.incontact];
np = length(incontact);
npic = sum(incontact);
dz0 = [contacts.dz];

%actuators
vis_act = u.vis_act; %abbreviate

%init output
out = ForwardDynOutput();
out.vc0 = NaN(3,np);
out.modelf_contact = zeros(3,np);

out.iter = 0;
out.cost = 0;

%END PRE-PROCESS


%compute tau
tau = zeros(nv,1);

%contact forces
if npic > 0
    if mdl.nw > 0
        Acontact = wheelJacobians(mdl, HT_world, contacts);
    elseif mdl.nt > 0
        Acontact = trackJacobians(mdl, HT_world, contacts);
    end
    
    vc0_incontact = reshape(Acontact*qvel0,3,npic);
    out.vc0(:,incontact) = vc0_incontact;

    Rw0 = radii(whichwt) .* qvel0(wsframeinds(whichwt)+5)';
    fc = feval(mdl.wgc_fh,mdl.wgc_p,incontact,vc0_incontact,Rw0(incontact),dz0(incontact));
    tau = tau + Acontact'*fc(:);
        
    out.modelf_contact(:,incontact) = fc;
end

%actuator forces
[out.modelf_act,u.err] = feval(mdl.act_fh, mdl.act_p, u.cmd, qvel0(vis_act), u.interr);
tau(vis_act) = tau(vis_act) + out.modelf_act;

%joint constraint forces
if ~isempty(mdl.hjc_fh)
    [~,Jc,out.modelf_joint] = feval(mdl.hjc_fh, mdl, state0(isjd), qvel0(7:end));
    tau(7:end) = tau(7:end) + Jc'*out.modelf_joint;
end

qacc = H\(tau-C);

