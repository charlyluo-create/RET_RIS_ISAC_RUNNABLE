function res = evaluate_fixed_ret_ris_config(cfg,scen,theta_deg,phi,gamma,P_total,use_ris)
%EVALUATE_FIXED_RET_RIS_CONFIG Evaluate a fixed RET/RIS configuration.
%
% RET and RIS remain fixed at theta_deg and phi, while the active
% communication and sensing beamformers are re-optimized from the current
% channel realization. This corresponds to the fast-timescale update in the
% multi-timescale control model.
% Optional USE_RIS defaults to true. Set false for RET-only/no-RIS ablation.

if nargin<7 || isempty(use_ris); use_ris=true; end

[~,ti]=min(abs(scen.theta_range-theta_deg));
[Hd_u,Hd_t,G]=apply_ret_channels(ti,scen.Hd_user,scen.Hd_tar,scen.G, ...
    scen.gain_u,scen.gain_t,scen.gain_ris);
res=evaluate_dynamic_isac(phi,use_ris,Hd_u,scen.Hr_user,Hd_t,scen.Hr_tar,G, ...
    gamma,sqrt(cfg.sigma_c2),cfg.sigma_r2,P_total,cfg.target_weights);
res.theta_deg=scen.theta_range(ti);
if use_ris; res.phi=normalize_phase(phi,scen.N); else; res.phi=[]; end
end

function res=evaluate_dynamic_isac(phi,use_ris,Hd_u,Hr_u,Hd_t,Hr_t,G, ...
    gamma,sigma_c,sigma_r2,P_total,omega)
res=empty_result();
if use_ris
    Hu=build_effective_user_channel(phi,true,Hd_u,Hr_u,G);
else
    Hu=Hd_u;
end
[Wc,p_comm,feas]=solve_comm_min_power(Hu,gamma,sigma_c);
if ~feas || ~isfinite(p_comm) || p_comm>P_total*(1+1e-6); return; end
C=build_radar_matrix(use_ris,phi,Hd_t,Hr_t,G,sigma_r2,omega);
[Wr,p_radar]=allocate_radar_power(Hu,C,P_total,p_comm);
W=[Wc,Wr];
res.feasible=true;
res.p_comm=p_comm;
res.p_radar=p_radar;
res.obj_snr=real(trace(W'*C*W));
res.Wc=Wc; res.Wr=Wr; res.Hu=Hu; res.C=C;
end

function [Hd_u,Hd_t,Gt]=apply_ret_channels(ti,Hd_u0,Hd_t0,G0,gain_u,gain_t,gain_ris)
Hd_u=Hd_u0; Hd_t=Hd_t0;
for k=1:size(Hd_u0,2)
    Hd_u(:,k)=sqrt(gain_u(ti,k))*Hd_u0(:,k);
end
for t=1:size(Hd_t0,2)
    Hd_t(:,t)=sqrt(gain_t(ti,t))*Hd_t0(:,t);
end
Gt=sqrt(gain_ris(ti))*G0;
end

function phi=normalize_phase(phi,N)
if isempty(phi)||numel(phi)~=N
    phi=ones(N,1);
else
    phi=phi(:);
end
bad=abs(phi)<1e-12; phi(bad)=1; phi=phi./abs(phi);
end

function res=empty_result()
res=struct('feasible',false,'p_comm',NaN,'p_radar',NaN,'obj_snr',NaN, ...
    'Wc',[],'Wr',[],'Hu',[],'C',[],'theta_deg',NaN,'phi',[]);
end
