function [res_best,history,meta] = optimize_ret_ris_block_coordinate( ...
    cfg,scen,gamma,P_total,theta_seed_deg,phi_seed)
%OPTIMIZE_RET_RIS_BLOCK_COORDINATE Fine same-model derivative-free baseline.
%
% The baseline starts from the same Stage-I point used by MA-RCG-MF in Exp.7,
% alternates a local RET search and a block-wise RIS phase-grid search, and
% evaluates every tested candidate with the exact communication-QoS solver
% and the exact sensing objective.  The zero phase offset is NOT re-evaluated
% because the current point is already retained as the incumbent.  The
% returned META structure records exact evaluations and completed sweeps so
% the runtime/complexity comparison can be audited directly.

N=scen.N;
phi=normalize_phase(phi_seed,N);
[~,theta_idx]=min(abs(scen.theta_range-theta_seed_deg));
meta=struct('exact_evals',0,'ret_evals',0,'ris_evals',0, ...
    'sweeps_completed',0,'accepted_updates',0,'blocks',NaN,'phase_levels',NaN);

res_best=evaluate_fixed_ret_ris_config(cfg,scen,scen.theta_range(theta_idx), ...
    phi,gamma,P_total);
meta.exact_evals=meta.exact_evals+1;
if ~res_best.feasible
    theta_idx=scen.ref_idx;
    phi=ones(N,1);
    res_best=evaluate_fixed_ret_ris_config(cfg,scen,scen.theta_range(theta_idx), ...
        phi,gamma,P_total);
    meta.exact_evals=meta.exact_evals+1;
end
history=struct('ObjSNR',[],'Pcomm',[],'ThetaDeg',[]);
if ~res_best.feasible; return; end
history=append_history(history,res_best,scen.theta_range(theta_idx));

num_blocks=max(1,min(N,round(cfg.exp7_bcs_blocks)));
levels=max(2,round(cfg.exp7_bcs_phase_levels));
num_sweeps=max(1,round(cfg.exp7_bcs_sweeps));
meta.blocks=num_blocks; meta.phase_levels=levels;
block_id=ceil((1:N)*num_blocks/N);
% Current point implicitly represents the zero offset; search only non-zero
% phase offsets so all timed exact evaluations do useful comparator work.
phase_offsets=2*pi*(1:levels-1)/levels;

for sweep=1:num_sweeps
    improved=false;

    % Common-RET local step.
    lo=max(1,theta_idx-cfg.tilt_local_half_window);
    hi=min(numel(scen.theta_range),theta_idx+cfg.tilt_local_half_window);
    cand_best=res_best; cand_idx=theta_idx;
    for ti=lo:hi
        if ti==theta_idx; continue; end
        rr=evaluate_fixed_ret_ris_config(cfg,scen,scen.theta_range(ti),phi, ...
            gamma,P_total);
        meta.exact_evals=meta.exact_evals+1; meta.ret_evals=meta.ret_evals+1;
        if rr.feasible && rr.obj_snr>cand_best.obj_snr*(1+1e-10)
            cand_best=rr; cand_idx=ti;
        end
    end
    if cand_best.obj_snr>res_best.obj_snr*(1+1e-10)
        res_best=cand_best; theta_idx=cand_idx; improved=true;
        meta.accepted_updates=meta.accepted_updates+1;
        history=append_history(history,res_best,scen.theta_range(theta_idx));
    end

    % RIS block-coordinate phase-offset step.
    for bb=1:num_blocks
        ids=find(block_id==bb);
        local_best=res_best; local_phi=phi;
        for qq=1:numel(phase_offsets)
            phi_q=phi;
            phi_q(ids)=phi_q(ids)*exp(1i*phase_offsets(qq));
            phi_q=normalize_phase(phi_q,N);
            rr=evaluate_fixed_ret_ris_config(cfg,scen, ...
                scen.theta_range(theta_idx),phi_q,gamma,P_total);
            meta.exact_evals=meta.exact_evals+1; meta.ris_evals=meta.ris_evals+1;
            if rr.feasible && rr.obj_snr>local_best.obj_snr*(1+1e-10)
                local_best=rr; local_phi=phi_q;
            end
        end
        if local_best.obj_snr>res_best.obj_snr*(1+1e-10)
            res_best=local_best; phi=local_phi; improved=true;
            meta.accepted_updates=meta.accepted_updates+1;
            history=append_history(history,res_best,scen.theta_range(theta_idx));
        end
    end

    meta.sweeps_completed=sweep;
    if ~improved; break; end
end
res_best.theta_deg=scen.theta_range(theta_idx);
res_best.phi=phi;
end

function h=append_history(h,r,theta)
h.ObjSNR(end+1)=r.obj_snr;
h.Pcomm(end+1)=r.p_comm;
h.ThetaDeg(end+1)=theta;
end

function phi=normalize_phase(phi,N)
if isempty(phi)||numel(phi)~=N; phi=ones(N,1); else; phi=phi(:); end
bad=abs(phi)<1e-12; phi(bad)=1; phi=phi./abs(phi);
end
