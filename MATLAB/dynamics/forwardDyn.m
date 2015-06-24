function [qacc, u, out] = forwardDyn(mdl, state0, qvel0, u, HT_parent, HT_world, contacts, dt)
%forward dynamics, compute joint space acceleration.
%use suffix 0 in variable names to distinguish values at time i from values at time i+1, as needed
%INPUTS
%mdl:       WmrModel object
%state0:    ns x 1, required only as input for joint constraint function (assumed to be consistent with HT_parent, HT_world)
%qvel0:     nv x 1, joint space velocity
%u:         ControllerIO object, use .cmd, .interr, .vis_act
%HT_parent: 4x4x nf, Homogeneous transforms to parent coords, HT_parent(:,:,i) = HT^(parent)_(frame i)
%HT_world:  4x4x nf, to world coords, HT_world(:,:,i) = HT^(world)_(frame i)
%contacts:  1 x nw array of WheelContactGeom objects, or
%           1 x nt array of TrackContactGeom objects
%dt:        integration step size (only needed if using constraints)
%OUTPUTS
%qacc:      ns x 1, joint space acceleration, d/dt qvel
%u:         ControllerIO object, .err updated
%out:       ForwardDynOutput object

if isempty(u.vis_act)
    nv = mdl.nf + 5;
    u.vis_act = false(1,nv);
    u.vis_act(mdl.actframeinds+5) = true;
end

[H, C] = HandC(mdl, HT_parent, qvel0);

if mdl.use_constraints
    %minimize the norm of H*qacc-(tau-C) or equivalently:
    %min 0.5*qacc'*H*qacc - (tau-C)'*qacc
    %subject to:
    %A*qacc = b
    %solve using Lagrange multipliers:
    %[H A'; A 0]*[qacc; lambda] = [tau-C; b];
    
    [A, cis] = constraintJacobians(mdl, state0, qvel0, u.vis_act, HT_world, contacts);
    
    if ~isempty(A)
        if isempty(mdl.wgc_fh)
            [qacc, u, out] = forwardDynErpCfm(mdl, state0, qvel0, u, contacts, dt, H, C, A, cis);
        else
            [qacc, u, out] = forwardDynForceBalance(mdl, state0, qvel0, u, contacts, dt, H, C, A, cis);
        end
        
        return
    end
    
end

%unconstrained
[qacc, u, out] = forwardDynUnc(mdl, state0, qvel0, u, HT_world, contacts, H, C);








