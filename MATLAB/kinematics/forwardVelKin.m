function  [qvel, vc] = forwardVelKin(mdl, state, u, vis_act, HT_world, contacts)
%INPUTS
%mdl:       WmrModel object
%state:     ns x 1 vector
%u:         na x 1 vector, actuated (or fixed) joint rates
%vis_act:   nv x 1 logical, which elements of qvel are actuated. sum(vis_act) = na. If empty [] use WmrModel
%HT_world:  4x4x nf matrix, HT_world(:,:,i) is transform from frame i to world coords. WmrModel and contact frames.
%contacts:  1 x nw array of WheelContactGeom objects, or
%           1 x nt array of TrackContactGeom objects
%OUTPUTS
%qvel:      nv x 1 vector, joint space velocity
%vc:        3 x np, velocity of contact pts with respect to world in contact coords, NaN if not in contact
%TODO, dqvel/du

%get from WmrModel object
nv = mdl.nf+5; %number of elements in qvel
[~,~,isjd] = stateIs(length(state)); %is joint displacement

if isempty(vis_act)
    vis_act = false(1,nv);
    vis_act(mdl.actframeinds+5) = true;
end

%contacts
incontact = [contacts.incontact];
np = length(incontact);
npic = sum(incontact); %number of points in contact

%check for holonomic joint constraints
if ~isempty(mdl.hjc_fh)
    [c, Jc] = feval(mdl.hjc_fh,mdl,state(isjd));
    njc = length(c); %number of (holonomic) joint constraints
    Ajoint = [zeros(njc,6), Jc]; %Jacobian wrt qvel
else
    njc = 0;
end
nc = 3*npic+njc; %number of constraints

%allocate for constraint equation, A*qvel = v
A = zeros(nc,nv); %including columns for fixed DOF
v = zeros(nc,1); %desired values, padded with NaN

%logicals, for indexing rows of A or v
cis.contact = false(nc,1); cis.contact(1:3*npic) = true; %is contact pt constraint
cis.vz = false(nc,1); cis.vz(3:3:3*npic) = true; %is contact pt vz constraint

if mdl.nw > 0
    Acontact = wheelJacobians(mdl, HT_world, contacts);
    vz = -([contacts.dz] - mdl.dz_target) ./ mdl.tc_z;
elseif mdl.nt > 0
    Acontact = trackJacobians(mdl, HT_world, contacts);
    vz = [];
    for tno = 1:mdl.nt
        vz = [vz, -(contacts(tno).dz - mdl.dz_target(tno)) / mdl.tc_z(tno)]; %#ok<AGROW>
    end
end
A(cis.contact,:) = Acontact;

%Baumgarte stabilization
v(cis.vz) = vz(incontact);

if njc > 0
    ijc = nc-njc+1; %row index of first holonomic joint constraint in A
    A(ijc:end,:) = Ajoint;
    v(ijc:end) = -c ./ mdl.tc_j'; %Baumgarte stabilization
end


%if wheel is not actuated and not in contact, make actuated to avoid rank deficiency
u_nv = NaN(nv,1); %padded
u_nv(vis_act) = u;

for wno = 1:mdl.nw
    vi = mdl.wheelframeinds(wno) + 5; %index in qvel
    if ~vis_act(vi) && ~incontact(wno)
        vis_act(vi) = true;
        u_nv(vi) = 0;
    end
end


%remove the cols of A corresponding to fixed rates to right hand side
%solve only for free DoFs
%Afree*qvel_free = b
Afree = A(:,~vis_act);
b = v - A(:,vis_act)*u;


qvel_free = Afree\b; %only the free elements

qvel = zeros(nv,1);
qvel(vis_act) = u;
qvel(~vis_act) = qvel_free;


if ~isempty(mdl.bsm_fh)
    %body slip model
    R_body_to_world = HT_world(1:3,1:3,1);

    vbody_slip = feval(mdl.bsm_fh,mdl.bsm_p,qvel(1:6),mdl.grav,mdl.Is,R_body_to_world);
    
    if 1
        %add directly to body frame velocity, faster but may cause holonomic constraint violation
        qvel(1:6) = qvel(1:6) + vbody_slip;
    else
        %instead of adding directly, convert to wheel slip, re-solve for qvel

        vc_slip = Acontact(:,1:6)*vbody_slip;
        vc_slip(3:3:end) = 0; %leave z components unchanged, holonomic constraints unaffected
        v_ = v;
        v_(cis.contact) = v_(cis.contact) + vc_slip;
        b = v_ - A(:,vis_fix)*qvel_fix;
        
        qvel_free = Afree\b;
        qvel(~vis_fix) = qvel_free;
    end
    
end

if nargout > 1
    %slip velocities for output
    vc = NaN(3,np);
    vc(:,incontact) = reshape(A(cis.contact,:)*qvel,[3 npic]);
end










