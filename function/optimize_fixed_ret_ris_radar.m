function [phi_best,res_best,refine_meta] = optimize_fixed_ret_ris_radar( ...
    cfg,scen,theta_deg,phi_seed,gamma,P_total)
%OPTIMIZE_FIXED_RET_RIS_RADAR Refresh RIS at fixed slow-timescale RET tilt.
%
% The RET angle is held fixed. RIS coefficients are optimized for the exact
% weighted radar SNR while communication QoS and total power constraints are
% enforced. Active beamformers are re-solved for every phase candidate.
% Experiment 6 retains the stale phase and may add deterministic current-
% snapshot candidates, while disabling unrelated random multi-start gains.

refine_meta=struct('screened',0,'exact_evals',0,'accepts',0);
[~,ti]=min(abs(scen.theta_range-theta_deg));
[Hd_u,Hd_t,G]=apply_ret_channels(ti,scen.Hd_user,scen.Hd_tar,scen.G, ...
    scen.gain_u,scen.gain_t,scen.gain_ris);

sigma_c=sqrt(cfg.sigma_c2);
sigma_r2=cfg.sigma_r2;
N=scen.N; K=scen.K; M=scen.M; T=scen.T;
num_trials=cfg.num_phase_trials;
if isfield(cfg,'exp6_fresh_ris_trials') && ~isempty(cfg.exp6_fresh_ris_trials)
    num_trials=cfg.exp6_fresh_ris_trials;
end
warm_only=isfield(cfg,'exp6_warm_start_only') && cfg.exp6_warm_start_only;
candidate_mode='warm_only';
if isfield(cfg,'exp6_fresh_candidate_mode') && ~isempty(cfg.exp6_fresh_candidate_mode)
    candidate_mode=lower(char(cfg.exp6_fresh_candidate_mode));
elseif ~warm_only
    candidate_mode='full';
end
allow_random=isfield(cfg,'exp6_fresh_allow_random') && cfg.exp6_fresh_allow_random;

[phi_best,res_best]=select_phase_radar(phi_seed,num_trials,Hd_u,scen.Hr_user, ...
    Hd_t,scen.Hr_tar,G,gamma,sigma_c,sigma_r2,P_total,cfg.target_weights, ...
    candidate_mode,allow_random);
if ~res_best.feasible
    res_best.theta_deg=scen.theta_range(ti);
    res_best.phi=phi_best;
    return;
end

rho=cfg.rho_init;
A=make_A(res_best.Hu,res_best.Wc);
for it=1:cfg.max_ao_iter
    old=res_best.obj_snr;
    try
        phi_new=call_opt_phi_radar_gain_warm(N,K,M,T,Hd_t,scen.Hr_tar, ...
            Hd_u,scen.Hr_user,G,res_best.Wc,A,rho,cfg.target_weights,phi_best);
        phi_new=normalize_phase(phi_new,N);
    catch ME
        warning('Fresh-RIS update stopped at iteration %d: %s',it,ME.message);
        break;
    end

    phase_cands={phi_new,normalize_phase(phi_best+phi_new,N)};
    cand=empty_result(); cand.obj_snr=-Inf; phi_accept=phi_best;
    for cc=1:numel(phase_cands)
        rr=evaluate_dynamic_isac(phase_cands{cc},true,Hd_u,scen.Hr_user, ...
            Hd_t,scen.Hr_tar,G,gamma,sigma_c,sigma_r2,P_total,cfg.target_weights);
        if rr.feasible && rr.obj_snr>cand.obj_snr
            cand=rr; phi_accept=phase_cands{cc};
        end
    end
    if ~cand.feasible || cand.obj_snr<=old*(1+1e-8); break; end

    phi_best=phi_accept; res_best=cand;
    A=opt_a_min_pc(K,M,A,res_best.Hu,res_best.Wc,rho,sigma_c,gamma);
    rho=rho/cfg.c_rho;
    aux_res=auxiliary_residual(A,res_best.Hu,res_best.Wc);
    rel=abs(res_best.obj_snr-old)/max(abs(old),1e-12);
    if rel<cfg.tol_obj && aux_res<cfg.tol_aux && it>2; break; end
end
if strcmp(candidate_mode,'full') && isfield(cfg,'fixed_ris_enable_exact_polish') && ...
        cfg.fixed_ris_enable_exact_polish && res_best.feasible
    [phi_best,res_best,refine_meta]=refine_fixed_ret_ris_exact( ...
        cfg,scen,scen.theta_range(ti),phi_best,res_best,gamma,P_total);
end
res_best.theta_deg=scen.theta_range(ti);
res_best.phi=phi_best;
end

function [phi_best,res_best]=select_phase_radar(phi_warm,num_trials, ...
    Hd_u,Hr_u,Hd_t,Hr_t,G,gamma,sigma_c,sigma_r2,P_total,omega, ...
    candidate_mode,allow_random)
N=size(G,1); T=size(Hd_t,2);
C=cell(0,1);
C{end+1}=normalize_phase(phi_warm,N);

% Deterministic candidates are physically tied to the current snapshot and
% do not introduce unrelated random-restart gains into the aging metric.
if contains(candidate_mode,'target') || strcmp(candidate_mode,'full')
    C{end+1}=target_focus_phase(Hd_t,Hr_t,G,omega,[]);
    center_t=ceil(T/2);
    C{end+1}=target_focus_phase(Hd_t,Hr_t,G,omega,center_t);
end
if contains(candidate_mode,'communication') || strcmp(candidate_mode,'full')
    C{end+1}=communication_focus_phase(Hd_u,Hr_u,G);
end
if strcmp(candidate_mode,'full')
    C{end+1}=ones(N,1);
    C{end+1}=exp(1i*2*pi*(0:N-1)'/N);
end
if allow_random
    while numel(C)<num_trials
        C{end+1}=exp(1i*2*pi*rand(N,1));
    end
end
% Remove exact duplicates to avoid repeated CVX solves.
C=unique_phase_candidates(C);
if num_trials>0 && numel(C)>num_trials
    C=C(1:num_trials);
end

phi_best=C{1}; res_best=empty_result(); res_best.obj_snr=-Inf;
for q=1:numel(C)
    rr=evaluate_dynamic_isac(C{q},true,Hd_u,Hr_u,Hd_t,Hr_t,G, ...
        gamma,sigma_c,sigma_r2,P_total,omega);
    if rr.feasible && rr.obj_snr>res_best.obj_snr
        res_best=rr; phi_best=C{q};
    end
end
end

function phi=target_focus_phase(Hd_t,Hr_t,G,omega,target_idx)
N=size(G,1); score=zeros(N,1); T=size(Hd_t,2);
if nargin<5||isempty(target_idx); ids=1:T; else; ids=target_idx; end
for t=ids
    gt=conj(G)*conj(Hd_t(:,t));
    f2=2*conj(Hr_t(:,t).*gt);
    if numel(ids)==1 && ~isempty(target_idx); wt=1; else; wt=omega(t); end
    score=score+wt*f2;
end
phi=phase_from_score(score,N);
end

function phi=communication_focus_phase(Hd_u,Hr_u,G)
N=size(G,1); score=zeros(N,1); K=size(Hd_u,2);
for k=1:K
    gk=conj(G)*conj(Hd_u(:,k));
    score=score+2*conj(Hr_u(:,k).*gk)/K;
end
phi=phase_from_score(score,N);
end

function phi=phase_from_score(score,N)
if isempty(score)||norm(score)<1e-14
    phi=ones(N,1);
else
    phi=exp(-1i*angle(score(:)));
end
phi=normalize_phase(phi,N);
end

function phi=call_opt_phi_radar_gain_warm(N,K,M,T,Hd_t,Hr_t,Hd_u,Hr_u,G,W,A,rho,omega,phi0)
try
    phi=opt_phi_radar_gain(N,K,M,T,Hd_t,Hr_t,Hd_u,Hr_u,G,W,A,rho,omega,phi0);
catch ME
    if contains(lower(ME.message),'too many input')
        phi=opt_phi_radar_gain(N,K,M,T,Hd_t,Hr_t,Hd_u,Hr_u,G,W,A,rho,omega);
    else
        rethrow(ME);
    end
end
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
if ~feas||~isfinite(p_comm)||p_comm>P_total*(1+1e-6); return; end
C=build_radar_matrix(use_ris,phi,Hd_t,Hr_t,G,sigma_r2,omega);
[Wr,p_radar]=allocate_radar_power(Hu,C,P_total,p_comm);
W=[Wc,Wr];
res.feasible=true; res.p_comm=p_comm; res.p_radar=p_radar;
res.obj_snr=real(trace(W'*C*W));
res.Wc=Wc; res.Wr=Wr; res.Hu=Hu; res.C=C;
end

function A=make_A(Hu,Wc)
K=size(Hu,2); A=zeros(K,K);
for k=1:K
    for j=1:K
        A(k,j)=Hu(:,k)'*Wc(:,j);
    end
end
end

function r=auxiliary_residual(A,Hu,Wc)
Aphys=make_A(Hu,Wc); d=A(:)-Aphys(:); r=max(abs(d));
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

function C=unique_phase_candidates(C)
keep=true(1,numel(C));
for ii=1:numel(C)
    if ~keep(ii); continue; end
    for jj=ii+1:numel(C)
        if keep(jj) && norm(C{ii}-C{jj})/sqrt(max(numel(C{ii}),1))<1e-10
            keep(jj)=false;
        end
    end
end
C=C(keep);
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
