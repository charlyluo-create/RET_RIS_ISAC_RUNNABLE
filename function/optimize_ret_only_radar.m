function [theta_best,res_best,history] = optimize_ret_only_radar(cfg,scen,gamma,P_total)
%OPTIMIZE_RET_ONLY_RADAR Exhaustive hardware-grid RET search for radar SNR.
% RIS is disabled. For every admissible common electronic downtilt, the
% minimum-power communication beamformer and null-space sensing beam are
% re-optimized. The exact weighted four-path radar objective selects the
% best feasible tilt. This is the same-objective RET-only ablation baseline.

res_best=struct('feasible',false,'p_comm',NaN,'p_radar',NaN, ...
    'obj_snr',NaN,'Wc',[],'Wr',[],'Hu',[],'C',[],'theta_deg',NaN,'phi',[]);
theta_best=NaN;
history.theta_deg=scen.theta_range(:);
history.obj_snr=nan(numel(scen.theta_range),1);
history.p_comm=nan(numel(scen.theta_range),1);
history.feasible=false(numel(scen.theta_range),1);

for ii=1:numel(scen.theta_range)
    theta=scen.theta_range(ii);
    rr=evaluate_fixed_ret_ris_config(cfg,scen,theta,[],gamma,P_total,false);
    history.feasible(ii)=rr.feasible;
    if rr.feasible
        history.obj_snr(ii)=rr.obj_snr;
        history.p_comm(ii)=rr.p_comm;
        if ~res_best.feasible || rr.obj_snr>res_best.obj_snr*(1+1e-12)
            res_best=rr;
            theta_best=theta;
        end
    end
end
end
