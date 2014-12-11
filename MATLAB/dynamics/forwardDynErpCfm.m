function [qacc, u, out] = forwardDynErpCfm( mdl, state0, qvel0, u, contacts, dt, inv_H, C, A, cis)
%Compute *constrained* forward dynamics using erp, cfm like Open Dynamics Engine
%TODO, bound contact forces

%PRE-PROCESS

%get from WmrModel object
ns = length(state0);
nv = mdl.nf+5; %number of elements in qvel
[~,~,isjd] = stateIs(ns);

%to compute wheel/sprocket radius * angular rate
if mdl.nw > 0
    whichwt = 1:mdl.nw; %which wheel/track does contact point belong to
elseif mdl.nt > 0
    whichwt = uncatIndexArray([contacts.np]);
end

%contacts
incontact = [contacts.incontact];
np = length(incontact);
npic = sum(incontact);
dz0 = [contacts.dz];

%actuators
vis_act = u.vis_act; %abbreviate
na = sum(vis_act);

%init output
u.err = zeros(na,1);

out = ForwardDynOutput();
out.vc0 = NaN(3,np);
out.modelf_contact = zeros(3,np);

out.iter = 0;
out.cost = 0;

%
nc = size(A,1);
ncc = sum(cis.contact);
nac = sum(cis.act);
njc = sum(cis.joint);

vc0_incontact = reshape(A(cis.contact,:)*qvel0, 3, npic);

%END PRE-PROCESS

tau = zeros(nv,1);

%set b
b = NaN(nc,1);
cfm_diag = zeros(nc,1);

if ncc > 0
    vc = zeros(3,np);
    vc(3,:) = -mdl.erp_z(whichwt) .* dz0/dt;
    vc_incontact = vc(:,incontact);
    temp = vc_incontact - vc0_incontact;

    b(cis.contact) = temp(:);

    temp = [mdl.fds_x(whichwt); mdl.fds_y(whichwt); mdl.cfm_z(whichwt)];
    temp = temp(:,incontact);
    cfm_diag(cis.contact) = temp(:);
end
if nac > 0
    b(cis.act) = u.cmd - qvel0(vis_act);
end
if njc > 0
    c = feval(mdl.hjc_fh, mdl, state0(isjd), qvel0(7:end)); %TODO, pass as input?
    b(cis.joint) = -mdl.erp_j .* c/dt - A(cis.joint,:)*qvel0;
    cfm_diag(cis.joint) = mdl.cfm_j;
end
b = b/dt;


%Hl*lambda = fl
Cfm = diag(cfm_diag/dt); %why divide by dt? because solving for acceleration (not velocity)
Hl = A*inv_H*A' + Cfm;
fl = b - A*inv_H*(tau-C);

%TODO, bound contact forces
%-normal force >= 0
%-lon,lat force <= lim (mu*normal force)

lambda = Hl\fl;

qacc = inv_H*(tau-C+A'*lambda);

%set ForwardDynOutput
if ncc > 0
    out.modelf_contact(:,incontact) = reshape(lambda(cis.contact),3,npic);
end

if nac > 0
    out.modelf_act = lambda(cis.act);
end

if njc > 0
    out.modelf_joint = lambda(cis.joint);
end

out.vc = NaN(3,np);
if ncc > 0
    %account for cfm
    vc_incontact = vc_incontact - reshape(cfm_diag(cis.contact).*lambda(cis.contact),3,npic);
    out.vc(:,incontact) = vc_incontact;
end


