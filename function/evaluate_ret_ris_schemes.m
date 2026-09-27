function out = evaluate_ret_ris_schemes(cfg,scen,gamma,P_total,opts)
%EVALUATE_RET_RIS_SCHEMES Evaluate the five unified RET-RIS ISAC schemes.
%
% Scheme order
%   1) Fixed RET, No RIS
%   2) Fixed RET, RIS-MinPc
%   3) RET-NoRIS-MinPc
%   4) RET-RIS-MinPc
%   5) MA-RCG-MF v2 (RCG Stage-I + sequential multi-fidelity ATR Stage-II)
%
% Optional fifth input opts may contain warm starts used by scalability
% experiments without changing the default behavior of Experiments 1--4/6:
%   phi_seed_fixed_minpc, phi_seed_minpc, phi_seed_radar,
%   theta_seed_minpc_deg,
%   theta_seed_radar_deg, num_phase_trials_override,
%   max_ao_iter_override, skip_scheme1, skip_preliminary_schemes,
%   preserve_rng_when_skipping.
%
% Besides Pcomm/Pradar/final radar SNR, this function reports:
%   G_eff  : maximum unit-power sensing gain in Null(Hu')
%   Jcomm  : sensing contribution of the communication waveforms
%   Jradar : sensing contribution of the dedicated radar waveform
% These quantities expose the two mechanisms investigated in the paper:
% communication-power release and sensing-channel shaping.  Scheme 5 also
% returns out.PreRefineProposed, i.e. the Stage-I RCG-AO solution
% immediately before the v2 sequential multi-fidelity / adaptive-trust-region
% exact-objective refinement.  This
% snapshot is obtained at zero extra solver cost and is used only in Exp.7
% to quantify the benefit of the proposed refinement stage.

if nargin < 5 || isempty(opts); opts=struct(); end

M=scen.M; N=scen.N; K=scen.K;
phi_seed_fixed_minpc=get_field_or(opts,'phi_seed_fixed_minpc',[]);
phi_seed_minpc=get_field_or(opts,'phi_seed_minpc',[]);
phi_seed_radar=get_field_or(opts,'phi_seed_radar',[]);
theta_seed_minpc_deg=get_field_or(opts,'theta_seed_minpc_deg',[]);
theta_seed_radar_deg=get_field_or(opts,'theta_seed_radar_deg',[]);
max_scheme=get_field_or(opts,'max_scheme',5);
max_scheme=max(1,min(5,round(max_scheme)));
num_phase_trials=max(4,round(get_field_or(opts,'num_phase_trials_override', ...
    cfg.num_phase_trials)));
max_ao_iter=max(1,round(get_field_or(opts,'max_ao_iter_override', ...
    cfg.max_ao_iter)));
skip_scheme1=logical(get_field_or(opts,'skip_scheme1',false));
skip_preliminary=logical(get_field_or(opts,'skip_preliminary_schemes',false));
preserve_rng_when_skipping=logical(get_field_or(opts, ...
    'preserve_rng_when_skipping',true));
if ~isempty(phi_seed_fixed_minpc); phi_seed_fixed_minpc=normalize_phase(phi_seed_fixed_minpc,N); end
if ~isempty(phi_seed_minpc); phi_seed_minpc=normalize_phase(phi_seed_minpc,N); end
if ~isempty(phi_seed_radar); phi_seed_radar=normalize_phase(phi_seed_radar,N); end
if skip_preliminary && (isempty(phi_seed_minpc) || isempty(theta_seed_minpc_deg))
    error(['skip_preliminary_schemes requires phi_seed_minpc and ' ...
        'theta_seed_minpc_deg.']);
end
theta_range=scen.theta_range;
numTheta=numel(theta_range);
sigma_c=sqrt(cfg.sigma_c2);
sigma_r2=cfg.sigma_r2;

out.Pcomm=nan(1,5); out.Pradar=nan(1,5); out.Obj=nan(1,5);
out.Feasible=false(1,5); out.Theta=nan(1,5);
out.G_eff=nan(1,5); out.Jcomm=nan(1,5); out.Jradar=nan(1,5);
out.Phi=cell(1,5);
out.History=cell(1,5);
out.RuntimeBlock=nan(1,5); out.RuntimeCumulative=nan(1,5);
out.SolverStatsBlock=cell(1,5); out.SolverStatsCumulative=cell(1,5);
out.PreRefineProposed=empty_snapshot();
out.ProposedMeta=struct();

% Common fixed-tilt channels.
[Hd_u_ref,Hd_t_ref,G_ref]=apply_ret_channels(scen.ref_idx, ...
    scen.Hd_user,scen.Hd_tar,scen.G, ...
    scen.gain_u,scen.gain_t,scen.gain_ris);

%% Scheme 1: fixed tilt, no RIS
t1=0;
if ~skip_scheme1
    ss0=comm_solver_telemetry('snapshot');
    tt=tic;
    r1=evaluate_dynamic_isac([],false,Hd_u_ref,scen.Hr_user, ...
        Hd_t_ref,scen.Hr_tar,G_ref,gamma,sigma_c,sigma_r2, ...
        P_total,cfg.target_weights);
    t1=toc(tt);
    out=store_result(out,1,r1,theta_range(scen.ref_idx),[]);
    out.History{1}=single_point_history(r1,theta_range(scen.ref_idx));
    ss1=comm_solver_telemetry('snapshot');
    out.SolverStatsBlock{1}=comm_solver_stats_delta(ss0,ss1);
    out.SolverStatsCumulative{1}=out.SolverStatsBlock{1};
else
    out.SolverStatsBlock{1}=zero_solver_stats();
    out.SolverStatsCumulative{1}=zero_solver_stats();
end

%% Scheme 2: fixed tilt, RIS-MinPc
t2=0; r2=empty_result(); phi2=phi_seed_fixed_minpc;
if skip_preliminary && preserve_rng_when_skipping
    % Scheme 2 creates max(1,num_phase_trials-4) random N-vectors during
    % candidate initialization. Consume the identical stream so Schemes 4--5
    % receive exactly the same random candidates as the full evaluator.
    rand(N,max(1,num_phase_trials-4)); %#ok<RAND>
end
if ~skip_preliminary
    ss0=comm_solver_telemetry('snapshot');
    tt=tic;
    if isempty(phi_seed_fixed_minpc); phi2_init=ones(N,1); else; phi2_init=phi_seed_fixed_minpc; end
    [phi2,comm2,h2]=optimize_phi_minpc(phi2_init,num_phase_trials, ...
        Hd_u_ref,scen.Hr_user,G_ref,gamma,sigma_c,P_total, ...
        max_ao_iter,cfg.rho_init,cfg.c_rho,cfg.tol_obj,cfg.tol_aux);
    r2=finalize_isac_from_comm(comm2,phi2,true,Hd_t_ref,scen.Hr_tar,G_ref, ...
        sigma_r2,P_total,cfg.target_weights);
    t2=toc(tt);
    out=store_result(out,2,r2,theta_range(scen.ref_idx),phi2);
    out.History{2}=h2;
    ss2=comm_solver_telemetry('snapshot');
    out.SolverStatsBlock{2}=comm_solver_stats_delta(ss0,ss2);
    out.SolverStatsCumulative{2}=out.SolverStatsBlock{2};
else
    out.SolverStatsBlock{2}=zero_solver_stats();
    out.SolverStatsCumulative{2}=zero_solver_stats();
end

%% Scheme 3: RET-NoRIS-MinPc
t3=0; comm3=empty_comm(); idx3=scen.ref_idx;
if skip_preliminary
    [~,idx3]=min(abs(theta_range-theta_seed_minpc_deg));
    out.SolverStatsBlock{3}=zero_solver_stats();
    out.SolverStatsCumulative{3}=zero_solver_stats();
else
    ss0=comm_solver_telemetry('snapshot');
    tt=tic;
    for ti=1:numTheta
        [Hd_u_t,~,~]=apply_ret_channels(ti,scen.Hd_user,scen.Hd_tar,scen.G, ...
            scen.gain_u,scen.gain_t,scen.gain_ris);
        cc=evaluate_comm([],false,Hd_u_t,scen.Hr_user,scen.G, ...
            gamma,sigma_c,P_total);
        if cc.feasible && (~comm3.feasible || cc.p_comm<comm3.p_comm)
            comm3=cc; idx3=ti;
        end
    end
    [~,Hd_t3,G3]=apply_ret_channels(idx3,scen.Hd_user,scen.Hd_tar,scen.G, ...
        scen.gain_u,scen.gain_t,scen.gain_ris);
    r3=finalize_isac_from_comm(comm3,[],false,Hd_t3,scen.Hr_tar,G3, ...
        sigma_r2,P_total,cfg.target_weights);
    t3=toc(tt);
    out=store_result(out,3,r3,theta_range(idx3),[]);
    out.History{3}=single_point_history(r3,theta_range(idx3));
    ss3=comm_solver_telemetry('snapshot');
    out.SolverStatsBlock{3}=comm_solver_stats_delta(ss0,ss3);
    out.SolverStatsCumulative{3}=out.SolverStatsBlock{3};
end

%% Scheme 4: RET-RIS-MinPc
ss0=comm_solver_telemetry('snapshot');
tt=tic;
if ~isempty(theta_seed_minpc_deg)
    [~,init_idx4]=min(abs(theta_range-theta_seed_minpc_deg));
elseif comm3.feasible
    init_idx4=idx3;
else
    init_idx4=scen.ref_idx;
end
if ~isempty(phi_seed_minpc)
    seed4=phi_seed_minpc;
elseif r2.feasible
    seed4=phi2;
else
    seed4=ones(N,1);
end
[phi4,idx4,comm4,h4]=optimize_joint_ret_ris_minpc(seed4,init_idx4, ...
    num_phase_trials,scen.Hd_user,scen.Hr_user,scen.G, ...
    scen.gain_u,scen.gain_ris,theta_range,gamma,sigma_c,P_total, ...
    max_ao_iter,cfg.tilt_local_half_window,cfg.rho_init, ...
    cfg.c_rho,cfg.tol_obj,cfg.tol_aux);
[~,Hd_t4,G4]=apply_ret_channels(idx4,scen.Hd_user,scen.Hd_tar,scen.G, ...
    scen.gain_u,scen.gain_t,scen.gain_ris);
r4=finalize_isac_from_comm(comm4,phi4,true,Hd_t4,scen.Hr_tar,G4, ...
    sigma_r2,P_total,cfg.target_weights);
t4=toc(tt);
out=store_result(out,4,r4,theta_range(idx4),phi4);
out.History{4}=h4;
ss4=comm_solver_telemetry('snapshot');
out.SolverStatsBlock{4}=comm_solver_stats_delta(ss0,ss4);
out.SolverStatsCumulative{4}=comm_solver_stats_add(out.SolverStatsBlock{2},out.SolverStatsBlock{3},out.SolverStatsBlock{4});

if max_scheme<5
    out.RuntimeBlock=[t1 t2 t3 t4 NaN];
    out.RuntimeCumulative=[t1, t2, t3, t2+t3+t4, NaN];
    out.SolverStatsBlock{5}=zero_solver_stats();
    out.SolverStatsCumulative{5}=out.SolverStatsCumulative{4};
    return;
end

%% Scheme 5: Proposed MA-RCG-MF RET-RIS-RadarSNR
ss0=comm_solver_telemetry('snapshot');
tt=tic;
if ~isempty(phi_seed_radar)
    seed5=phi_seed_radar;
elseif r4.feasible
    seed5=phi4;
elseif r2.feasible
    seed5=phi2;
else
    seed5=ones(N,1);
end
if ~isempty(theta_seed_radar_deg)
    [~,radar_primary_idx]=min(abs(theta_range-theta_seed_radar_deg));
else
    radar_primary_idx=scen.target_idx;
end
prop_opts=build_proposed_options(cfg);
[phi5,idx5,r5,h5,meta5]=optimize_joint_ret_ris_radar(seed5,radar_primary_idx, ...
    [idx4 scen.ref_idx],num_phase_trials,scen.Hd_user, ...
    scen.Hr_user,scen.Hd_tar,scen.Hr_tar,scen.G, ...
    scen.gain_u,scen.gain_t,scen.gain_ris,theta_range, ...
    gamma,sigma_c,sigma_r2,P_total,cfg.target_weights, ...
    max_ao_iter,cfg.tilt_local_half_window,cfg.rho_init, ...
    cfg.c_rho,cfg.tol_obj,cfg.tol_aux,prop_opts);
t5=toc(tt);
out=store_result(out,5,r5,theta_range(idx5),phi5);
out.History{5}=h5;
out.ProposedMeta=meta5;
ss5=comm_solver_telemetry('snapshot');
out.SolverStatsBlock{5}=comm_solver_stats_delta(ss0,ss5);
out.SolverStatsCumulative{5}=comm_solver_stats_add(out.SolverStatsBlock{2},out.SolverStatsBlock{3},out.SolverStatsBlock{4},out.SolverStatsBlock{5});
if isfield(meta5,'stage1_result') && meta5.stage1_result.feasible
    out.PreRefineProposed=make_snapshot(meta5.stage1_result, ...
        theta_range(meta5.stage1_idx),meta5.stage1_phi,meta5.stage1_history, ...
        meta5.stage1_runtime_s,t2+t3+t4+meta5.stage1_runtime_s, ...
        comm_solver_stats_add(out.SolverStatsCumulative{4},meta5.stage1_solver_stats));
end

out.RuntimeBlock=[t1 t2 t3 t4 t5];
% Proposed joint schemes use the fixed-RIS and RET-only solutions as
% physically meaningful initializers. The cumulative time reflects the
% complete executable path of each scheme.
out.RuntimeCumulative=[t1, t2, t3, t2+t3+t4, t2+t3+t4+t5];
end

function S=zero_solver_stats()
S=comm_solver_telemetry('snapshot'); fn=fieldnames(S);
for ii=1:numel(fn); if isnumeric(S.(fn{ii})); S.(fn{ii})=0; end; end
end

function value=get_field_or(s,name,default)
if isstruct(s) && isfield(s,name) && ~isempty(s.(name))
    value=s.(name);
else
    value=default;
end
end

function opts=build_proposed_options(cfg)
opts=struct();
opts.mechanism_fallback=get_cfg_or(cfg,'proposed_mechanism_fallback',true);
opts.hybrid_comm_weight_min=get_cfg_or(cfg,'proposed_hybrid_comm_weight_min',0.2);
opts.hybrid_comm_weight_max=get_cfg_or(cfg,'proposed_hybrid_comm_weight_max',0.8);
% Stage-I multi-fidelity screening: RCG generates several geodesic candidates,
% the mechanism proxy ranks all RET/RIS combinations, and only a small adaptive
% top-K subset is passed to the exact inner communication solver.
opts.stage1_geodesic_fractions=get_cfg_or(cfg,'proposed_stage1_geodesic_fractions',[1 0.5]);
opts.stage1_exact_topk_initial=get_cfg_or(cfg,'proposed_stage1_exact_topk_initial',2);
opts.stage1_exact_topk_max=get_cfg_or(cfg,'proposed_stage1_exact_topk_max',4);
opts.stage1_expand_on_stall=get_cfg_or(cfg,'proposed_stage1_expand_on_stall',true);
% Stage-II adaptive coarse-to-fine active-block refinement.
opts.polish_blocks_schedule=get_cfg_or(cfg,'proposed_polish_blocks_schedule',[8 16]);
opts.polish_phase_step_schedule=get_cfg_or(cfg,'proposed_polish_phase_step_schedule', ...
    {[22.5 11.25],[11.25 5.625]});
opts.polish_active_blocks_schedule=get_cfg_or(cfg,'proposed_polish_active_blocks_schedule', ...
    opts.polish_blocks_schedule);
opts.polish_exact_topk_initial_schedule=get_cfg_or(cfg,'proposed_polish_exact_topk_initial_schedule',[2 2]);
opts.polish_exact_topk_max_schedule=get_cfg_or(cfg,'proposed_polish_exact_topk_max_schedule',[4 4]);
opts.polish_ret_half_window=get_cfg_or(cfg,'proposed_polish_ret_half_window',1);
opts.polish_ret_exact_topk_initial=get_cfg_or(cfg,'proposed_polish_ret_exact_topk_initial',1);
opts.polish_ret_exact_topk_max=get_cfg_or(cfg,'proposed_polish_ret_exact_topk_max',2);
opts.polish_try_fine_on_stall=get_cfg_or(cfg,'proposed_polish_try_fine_on_stall',true);
opts.polish_min_gain_dB=get_cfg_or(cfg,'proposed_polish_min_gain_dB',0);
opts.refine_accept_tol=get_cfg_or(cfg,'proposed_refine_accept_tol',1e-10);

% MA-RCG-MF v2 Stage-II options.  The R3 diagnostic showed that the main
% loss was candidate-generation capacity, not Top-K screening.  V2 therefore
% visits every block sequentially and adapts an individual trust radius.
opts.stage2_variant=get_cfg_or(cfg,'proposed_stage2_variant','sequential_mf_atr');
opts.v2_blocks=get_cfg_or(cfg,'v2_blocks',8);
opts.v2_max_sweeps=get_cfg_or(cfg,'v2_max_sweeps',3);
opts.v2_global_probe_levels=get_cfg_or(cfg,'v2_global_probe_levels',16);
opts.v2_global_probe_sweeps=get_cfg_or(cfg,'v2_global_probe_sweeps',1);
opts.v2_trust_init_deg=get_cfg_or(cfg,'v2_trust_init_deg',45);
opts.v2_trust_min_deg=get_cfg_or(cfg,'v2_trust_min_deg',2.8125);
opts.v2_trust_max_deg=get_cfg_or(cfg,'v2_trust_max_deg',135);
opts.v2_trust_expand=get_cfg_or(cfg,'v2_trust_expand',1.5);
opts.v2_trust_shrink=get_cfg_or(cfg,'v2_trust_shrink',0.5);
opts.v2_local_fractions=get_cfg_or(cfg,'v2_local_fractions',[1 0.5 0.25]);
opts.v2_exact_topk_initial=get_cfg_or(cfg,'v2_exact_topk_initial',2);
opts.v2_exact_topk_max=get_cfg_or(cfg,'v2_exact_topk_max',4);
opts.v2_expand_budget_on_stall=get_cfg_or(cfg,'v2_expand_budget_on_stall',true);
opts.v2_proxy_tie_gap_dB=get_cfg_or(cfg,'v2_proxy_tie_gap_dB',0.05);
opts.v2_ret_half_window=get_cfg_or(cfg,'v2_ret_half_window',1);
opts.v2_ret_exact_topk_initial=get_cfg_or(cfg,'v2_ret_exact_topk_initial',1);
opts.v2_ret_exact_topk_max=get_cfg_or(cfg,'v2_ret_exact_topk_max',2);
opts.v2_min_gain_dB=get_cfg_or(cfg,'v2_min_gain_dB',2e-3);
opts.v2_sweep_stop_gain_dB=get_cfg_or(cfg,'v2_sweep_stop_gain_dB',5e-3);
opts.v2_max_noaccept_sweeps=get_cfg_or(cfg,'v2_max_noaccept_sweeps',2);
end

function value=get_cfg_or(cfg,name,default)
if isstruct(cfg) && isfield(cfg,name) && ~isempty(cfg.(name)); value=cfg.(name); else; value=default; end
end

function meta=empty_proposed_meta()
meta=struct('stage1_runtime_s',NaN,'refinement_runtime_s',0, ...
    'stage1_result',empty_result(),'stage1_phi',[],'stage1_idx',NaN, ...
    'stage1_history',empty_history(),'stage1_obj_snr',NaN,'stage1_p_comm',NaN, ...
    'stage1_screened',0,'stage1_exact_evals',0,'stage1_accepts',0, ...
    'stage1_solver_stats',zero_solver_stats(),'refinement_solver_stats',zero_solver_stats(), ...
    'mechanism_evals',0,'mechanism_accepted',false,'polish_evals',0, ...
    'polish_accepts',0,'polish_screened',0,'polish_round_accepts',[], ...
    'stage2_variant','','v2_sweeps_completed',0,'v2_block_visits',0, ...
    'v2_budget_expansions',0,'v2_ret_evals',0,'v2_ris_evals',0, ...
    'v2_sweep_gain_dB',[],'v2_sweep_accepts',[],'v2_trust_final_deg',[], ...
    'v2_block_accept_count',[],'v2_block_exact_evals',[], ...
    'total_exact_refine_evals',0,'total_exact_algorithm_evals',0,'final_obj_snr',NaN, ...
    'final_p_comm',NaN,'refinement_gain_dB',NaN,'monotonicity_guard_triggered',false);
end

function s=empty_snapshot()
s=struct('Pcomm',NaN,'Pradar',NaN,'Obj',NaN,'Feasible',false, ...
    'Theta',NaN,'Phi',[],'G_eff',NaN,'Jcomm',NaN,'Jradar',NaN, ...
    'History',empty_history(),'RuntimeBlock',NaN,'RuntimeCumulative',NaN, ...
    'SolverStatsCumulative',zero_solver_stats());
end

function s=make_snapshot(res,theta,phi,history,tblock,tcumulative,solverStats)
s=empty_snapshot();
if ~res.feasible; return; end
s.Pcomm=res.p_comm; s.Pradar=res.p_radar; s.Obj=res.obj_snr;
s.Feasible=true; s.Theta=theta; s.Phi=phi;
s.G_eff=effective_sensing_gain(res.Hu,res.C);
s.Jcomm=real(trace(res.Wc'*res.C*res.Wc));
s.Jradar=real(res.Wr'*res.C*res.Wr);
s.History=history; s.RuntimeBlock=tblock; s.RuntimeCumulative=tcumulative;
if nargin>=7 && ~isempty(solverStats); s.SolverStatsCumulative=solverStats; end
end

function out=store_result(out,s,res,theta,phi)
if ~res.feasible; return; end
out.Pcomm(s)=res.p_comm;
out.Pradar(s)=res.p_radar;
out.Obj(s)=res.obj_snr;
out.Feasible(s)=true;
out.Theta(s)=theta;
out.Phi{s}=phi;
out.G_eff(s)=effective_sensing_gain(res.Hu,res.C);
out.Jcomm(s)=real(trace(res.Wc'*res.C*res.Wc));
out.Jradar(s)=real(res.Wr'*res.C*res.Wr);
end

function g=effective_sensing_gain(Hu,C)
% Unit-radar-power gain under the communication-null-space constraint.
Z=null(Hu');
if isempty(Z)
    g=0;
    return;
end
B=(Z'*C*Z); B=(B+B')/2;
ev=real(eig(B));
g=max(max(ev),0);
end

%% ======================== Min-Pc Optimization =============================
function [phi_best,comm_best,history]=optimize_phi_minpc(phi_warm,num_trials, ...
    Hd_u,Hr_u,G,gamma,sigma_c,P_total,max_iter,rho_init,c_rho,tol_obj,tol_aux)

N=size(G,1); K=size(Hd_u,2); M=size(Hd_u,1);
[phi_best,comm_best]=select_phase_minpc(phi_warm,num_trials,Hd_u,Hr_u,G, ...
    gamma,sigma_c,P_total);
history=empty_history();
if comm_best.feasible
    history=append_history(history,comm_best.p_comm,NaN,NaN,NaN);
end
if ~comm_best.feasible; return; end
rho=rho_init; A=make_A(comm_best.Hu,comm_best.Wc);

for it=1:max_iter
    old=comm_best.p_comm;
    try
        phi_new=call_opt_phi_minpc_warm(N,K,M,Hd_u,Hr_u,G, ...
            comm_best.Wc,A,rho,phi_best);
        phi_new=normalize_phase(phi_new,N);
    catch
        break;
    end

    phase_cands={phi_new,normalize_phase(phi_best+phi_new,N)};
    cand=empty_comm(); cand_phi=phi_best;
    for cc=1:numel(phase_cands)
        rr=evaluate_comm(phase_cands{cc},true,Hd_u,Hr_u,G,gamma,sigma_c,P_total);
        if rr.feasible && (~cand.feasible || rr.p_comm<cand.p_comm)
            cand=rr; cand_phi=phase_cands{cc};
        end
    end
    if ~cand.feasible || cand.p_comm>=old*(1-1e-8); break; end

    phi_best=cand_phi; comm_best=cand;
    A=opt_a_min_pc(K,M,A,comm_best.Hu,comm_best.Wc,rho,sigma_c,gamma);
    rho=rho/c_rho;
    aux_res=auxiliary_residual(A,comm_best.Hu,comm_best.Wc);
    rel=abs(comm_best.p_comm-old)/max(old,1e-12);
    history=append_history(history,comm_best.p_comm,NaN,NaN,aux_res);
    if rel<tol_obj && aux_res<tol_aux && it>2; break; end
end
end

function [phi_best,best_idx,comm_best,history]=optimize_joint_ret_ris_minpc( ...
    phi_seed,init_idx,num_trials,Hd_u0,Hr_u,G0,gain_u,gain_ris,theta_range, ...
    gamma,sigma_c,P_total,max_iter,half_window,rho_init,c_rho,tol_obj,tol_aux)

[Hd_u,G]=apply_ret_comm_channels(init_idx,Hd_u0,G0,gain_u,gain_ris);
[phi_best,comm_best]=select_phase_minpc(phi_seed,num_trials,Hd_u,Hr_u,G, ...
    gamma,sigma_c,P_total);
best_idx=init_idx;
history=empty_history();
if comm_best.feasible
    history=append_history(history,comm_best.p_comm,NaN,theta_range(best_idx),NaN);
end
if ~comm_best.feasible; return; end
N=size(G,1); K=size(Hd_u,2); M=size(Hd_u,1);
rho=rho_init; A=make_A(comm_best.Hu,comm_best.Wc);

for it=1:max_iter
    old=comm_best.p_comm;
    try
        phi_new=call_opt_phi_minpc_warm(N,K,M,Hd_u,Hr_u,G, ...
            comm_best.Wc,A,rho,phi_best);
        phi_new=normalize_phase(phi_new,N);
    catch
        break;
    end

    lo=max(1,best_idx-half_window); hi=min(numel(theta_range),best_idx+half_window);
    phase_cands={phi_new,normalize_phase(phi_best+phi_new,N)};
    cand=empty_comm(); cand_idx=best_idx; cand_phi=phi_best;
    cand_Hd_u=Hd_u; cand_G=G;
    for cc=1:numel(phase_cands)
        for ti=lo:hi
            [Hd_u_t,G_t]=apply_ret_comm_channels(ti,Hd_u0,G0,gain_u,gain_ris);
            rr=evaluate_comm(phase_cands{cc},true,Hd_u_t,Hr_u,G_t, ...
                gamma,sigma_c,P_total);
            if rr.feasible && (~cand.feasible || rr.p_comm<cand.p_comm)
                cand=rr; cand_idx=ti; cand_phi=phase_cands{cc};
                cand_Hd_u=Hd_u_t; cand_G=G_t;
            end
        end
    end
    if ~cand.feasible || cand.p_comm>=old*(1-1e-8); break; end

    phi_best=cand_phi; best_idx=cand_idx; comm_best=cand;
    Hd_u=cand_Hd_u; G=cand_G;
    A=opt_a_min_pc(K,M,A,comm_best.Hu,comm_best.Wc,rho,sigma_c,gamma);
    rho=rho/c_rho;
    aux_res=auxiliary_residual(A,comm_best.Hu,comm_best.Wc);
    rel=abs(comm_best.p_comm-old)/max(old,1e-12);
    history=append_history(history,comm_best.p_comm,NaN,theta_range(best_idx),aux_res);
    if rel<tol_obj && aux_res<tol_aux && it>2; break; end
end
end

function [phi_best,comm_best]=select_phase_minpc(phi_warm,num_trials, ...
    Hd_u,Hr_u,G,gamma,sigma_c,P_total)
N=size(G,1);
C=cell(0,1);
C{end+1}=normalize_phase(phi_warm,N);
C{end+1}=communication_focus_phase(Hd_u,Hr_u,G);
C{end+1}=ones(N,1);
C{end+1}=exp(1i*2*pi*(0:N-1)'/N);
C{end+1}=exp(1i*2*pi*rand(N,1));
while numel(C)<num_trials; C{end+1}=exp(1i*2*pi*rand(N,1)); end
C=C(1:min(numel(C),max(num_trials,4)));

phi_best=C{1}; comm_best=empty_comm();
for q=1:numel(C)
    rr=evaluate_comm(C{q},true,Hd_u,Hr_u,G,gamma,sigma_c,P_total);
    if rr.feasible && (~comm_best.feasible || rr.p_comm<comm_best.p_comm)
        comm_best=rr; phi_best=C{q};
    end
end
end

%% ======================== Radar-SNR Optimization ==========================
function [phi_best,best_idx,res_best,history,meta]=optimize_joint_ret_ris_radar( ...
    phi_seed,primary_idx,fallback_idx,num_trials,Hd_u0,Hr_u,Hd_t0,Hr_t,G0, ...
    gain_u,gain_t,gain_ris,theta_range,gamma,sigma_c,sigma_r2,P_total,omega, ...
    max_iter,half_window,rho_init,c_rho,tol_obj,tol_aux,prop_opts)
%OPTIMIZE_JOINT_RET_RIS_RADAR Fast mechanism-aware RCG/AO optimizer.
%
% The algorithm follows the performance/complexity principles used in recent
% high-level RIS-ISAC optimization work (AO + manifold/RCG + surrogate/MM),
% but adds two system-specific accelerations:
%   (i) adaptive multi-fidelity exact evaluation in Stage I; and
%   (ii) sequential multi-fidelity block refinement with adaptive trust regions in Stage II.
% All updates are finally accepted using the true structured radar objective.

if nargin<25 || isempty(prop_opts); prop_opts=struct(); end
stage1_timer=tic;
stage1_solver_start=comm_solver_telemetry('snapshot');
idx_cands=unique([primary_idx,fallback_idx(:).'],'stable');
res_best=empty_result(); res_best.obj_snr=-Inf;
phi_best=normalize_phase(phi_seed,size(G0,1)); best_idx=idx_cands(1);
Hd_u=[]; Hd_t=[]; G=[];
for qq=1:numel(idx_cands)
    ti=idx_cands(qq);
    [Hu0,Ht0,Gt0]=apply_ret_channels(ti,Hd_u0,Hd_t0,G0,gain_u,gain_t,gain_ris);
    [phi0,r0]=select_phase_radar(phi_seed,num_trials,Hu0,Hr_u,Ht0,Hr_t,Gt0, ...
        gamma,sigma_c,sigma_r2,P_total,omega);
    if r0.feasible && r0.obj_snr>res_best.obj_snr
        res_best=r0; phi_best=phi0; best_idx=ti; Hd_u=Hu0; Hd_t=Ht0; G=Gt0;
    end
end
history=empty_history(); meta=empty_proposed_meta();
if res_best.feasible
    history=append_history(history,res_best.p_comm,res_best.obj_snr,theta_range(best_idx),NaN);
end
if ~res_best.feasible
    meta.stage1_runtime_s=toc(stage1_timer);
    meta.stage1_result=res_best; meta.stage1_phi=phi_best; meta.stage1_idx=best_idx;
    meta.stage1_history=history;
    meta.stage1_solver_stats=comm_solver_stats_delta(stage1_solver_start,comm_solver_telemetry('snapshot'));
    return;
end

N=size(G,1); K=size(Hd_u,2); M=size(Hd_u,1); T=size(Hd_t,2);
rho=rho_init; A=make_A(res_best.Hu,res_best.Wc);
fractions=get_field_or(prop_opts,'stage1_geodesic_fractions',[1 0.5]);
fractions=unique(min(max(real(fractions(:).'),0.05),1),'stable');
q0=max(1,round(get_field_or(prop_opts,'stage1_exact_topk_initial',2)));
qmax=max(q0,round(get_field_or(prop_opts,'stage1_exact_topk_max',4)));
expand_on_stall=logical(get_field_or(prop_opts,'stage1_expand_on_stall',true));

% -------------------------- Stage I --------------------------------------
% RCG generates a high-quality unit-modulus search point. Several points on
% the geodesic from the current phase to the RCG point are combined with the
% neighboring RET grid. The cheap mechanism proxy ranks all combinations.
% Only adaptive top-K candidates invoke the exact communication solver.
for it=1:max_iter
    old=res_best.obj_snr;
    try
        phi_rcg=call_opt_phi_radar_gain_warm(N,K,M,T,Hd_t,Hr_t,Hd_u,Hr_u,G, ...
            res_best.Wc,A,rho,omega,phi_best);
        phi_rcg=normalize_phase(phi_rcg,N);
    catch
        break;
    end

    phase_cands=cell(1,numel(fractions));
    for ff=1:numel(fractions)
        phase_cands{ff}=phase_geodesic_blend(phi_best,phi_rcg,fractions(ff));
    end
    lo=max(1,best_idx-half_window); hi=min(numel(theta_range),best_idx+half_window);
    cand_phi=cell(0,1); cand_idx=[]; cand_score=[];
    for cc=1:numel(phase_cands)
        for ti=lo:hi
            [Hu_t,Ht_t,G_t]=apply_ret_channels(ti,Hd_u0,Hd_t0,G0,gain_u,gain_t,gain_ris);
            sc=mechanism_proxy(phase_cands{cc},Hu_t,Hr_u,Ht_t,Hr_t,G_t,res_best, ...
                gamma,sigma_c,sigma_r2,P_total,omega,prop_opts);
            cand_phi{end+1}=phase_cands{cc}; %#ok<AGROW>
            cand_idx(end+1)=ti; %#ok<AGROW>
            cand_score(end+1)=sc; %#ok<AGROW>
            meta.stage1_screened=meta.stage1_screened+1;
        end
    end
    [~,ord]=sort(cand_score,'descend');
    [cand,cand_phi_best,cand_idx_best,used,improved]=evaluate_ranked_joint_candidates( ...
        ord,cand_score,cand_phi,cand_idx,q0,res_best,Hd_u0,Hr_u,Hd_t0,Hr_t,G0, ...
        gain_u,gain_t,gain_ris,gamma,sigma_c,sigma_r2,P_total,omega,prop_opts);
    meta.stage1_exact_evals=meta.stage1_exact_evals+used;
    if ~improved && expand_on_stall && qmax>q0
        [cand2,phi2,idx2,used2,improved2]=evaluate_ranked_joint_candidates( ...
            ord,cand_score,cand_phi,cand_idx,qmax,res_best,Hd_u0,Hr_u,Hd_t0,Hr_t,G0, ...
            gain_u,gain_t,gain_ris,gamma,sigma_c,sigma_r2,P_total,omega,prop_opts,q0);
        meta.stage1_exact_evals=meta.stage1_exact_evals+used2;
        if improved2; cand=cand2; cand_phi_best=phi2; cand_idx_best=idx2; improved=true; end
    end
    if ~improved; break; end

    phi_best=cand_phi_best; best_idx=cand_idx_best; res_best=cand;
    [Hd_u,Hd_t,G]=apply_ret_channels(best_idx,Hd_u0,Hd_t0,G0,gain_u,gain_t,gain_ris);
    A=opt_a_min_pc(K,M,A,res_best.Hu,res_best.Wc,rho,sigma_c,gamma);
    rho=rho/c_rho;
    aux_res=auxiliary_residual(A,res_best.Hu,res_best.Wc);
    rel=abs(res_best.obj_snr-old)/max(abs(old),1e-12);
    history=append_history(history,res_best.p_comm,res_best.obj_snr,theta_range(best_idx),aux_res);
    meta.stage1_accepts=meta.stage1_accepts+1;
    if rel<tol_obj && aux_res<tol_aux && it>2; break; end
end

meta.stage1_runtime_s=toc(stage1_timer);
meta.stage1_result=res_best; meta.stage1_phi=phi_best; meta.stage1_idx=best_idx;
meta.stage1_history=history; meta.stage1_obj_snr=res_best.obj_snr;
meta.stage1_p_comm=res_best.p_comm;
meta.stage1_solver_stats=comm_solver_stats_delta(stage1_solver_start,comm_solver_telemetry('snapshot'));

% -------------------------- Stage II-A -----------------------------------
% One physically structured communication/radar geodesic escape candidate.
refine_timer=tic;
refine_solver_start=comm_solver_telemetry('snapshot');
if get_field_or(prop_opts,'mechanism_fallback',true)
    phi_tar=target_focus_phase(Hd_t,Hr_t,G,omega,[]);
    phi_com=communication_focus_phase(Hd_u,Hr_u,G);
    alpha=res_best.p_comm/max(P_total,eps);
    alpha=min(max(alpha,get_field_or(prop_opts,'hybrid_comm_weight_min',0.2)), ...
        get_field_or(prop_opts,'hybrid_comm_weight_max',0.8));
    phi_mech=phase_geodesic_blend(phi_tar,phi_com,alpha);
    rr=evaluate_dynamic_isac(phi_mech,true,Hd_u,Hr_u,Hd_t,Hr_t,G, ...
        gamma,sigma_c,sigma_r2,P_total,omega);
    meta.mechanism_evals=meta.mechanism_evals+1;
    if exact_improvement(rr,res_best,prop_opts)
        phi_best=phi_mech; res_best=rr; A=make_A(res_best.Hu,res_best.Wc);
        history=append_history(history,res_best.p_comm,res_best.obj_snr, ...
            theta_range(best_idx),auxiliary_residual(A,res_best.Hu,res_best.Wc));
        meta.mechanism_accepted=true;
    end
end

% -------------------------- Stage II-B -----------------------------------
% MA-RCG-MF v2: Sequential Multi-Fidelity Block Refinement + Adaptive
% Trust Region (SMFBR-ATR). Every RIS block is visited, so the active-block
% exclusion failure observed in the R3 diagnostic is removed.  Each accepted
% block update immediately becomes the incumbent for the next block, enabling
% cooperative multi-block gains.  A full-circle proxy-only probe is used in
% the first sweep; later sweeps use per-block local trust regions.  Only a
% small adaptive Top-K receives the exact communication/QoS solve.
stage2_variant=lower(char(get_field_or(prop_opts,'stage2_variant','sequential_mf_atr')));
if strcmp(stage2_variant,'sequential_mf_atr')
    [phi_best,best_idx,res_best,history,polish_meta]= ...
        sequential_multifidelity_trust_region_polish( ...
        phi_best,best_idx,res_best,history,Hd_u0,Hr_u,Hd_t0,Hr_t,G0, ...
        gain_u,gain_t,gain_ris,theta_range,gamma,sigma_c,sigma_r2,P_total,omega,prop_opts);
else
    [phi_best,best_idx,res_best,history,polish_meta]=adaptive_multifidelity_polish( ...
        phi_best,best_idx,res_best,history,Hd_u0,Hr_u,Hd_t0,Hr_t,G0, ...
        gain_u,gain_t,gain_ris,theta_range,gamma,sigma_c,sigma_r2,P_total,omega,prop_opts);
end
meta.stage2_variant=stage2_variant;
meta.polish_evals=polish_meta.evals; meta.polish_accepts=polish_meta.accepts;
meta.polish_screened=polish_meta.screened;
if isfield(polish_meta,'round_accepts'); meta.polish_round_accepts=polish_meta.round_accepts; end
if isfield(polish_meta,'sweeps_completed'); meta.v2_sweeps_completed=polish_meta.sweeps_completed; end
if isfield(polish_meta,'block_visits'); meta.v2_block_visits=polish_meta.block_visits; end
if isfield(polish_meta,'budget_expansions'); meta.v2_budget_expansions=polish_meta.budget_expansions; end
if isfield(polish_meta,'ret_evals'); meta.v2_ret_evals=polish_meta.ret_evals; end
if isfield(polish_meta,'ris_evals'); meta.v2_ris_evals=polish_meta.ris_evals; end
if isfield(polish_meta,'sweep_gain_dB'); meta.v2_sweep_gain_dB=polish_meta.sweep_gain_dB; end
if isfield(polish_meta,'sweep_accepts'); meta.v2_sweep_accepts=polish_meta.sweep_accepts; end
if isfield(polish_meta,'trust_final_deg'); meta.v2_trust_final_deg=polish_meta.trust_final_deg; end
if isfield(polish_meta,'block_accept_count'); meta.v2_block_accept_count=polish_meta.block_accept_count; end
if isfield(polish_meta,'block_exact_evals'); meta.v2_block_exact_evals=polish_meta.block_exact_evals; end

if meta.stage1_result.feasible && (~res_best.feasible || ...
        res_best.obj_snr<meta.stage1_result.obj_snr*(1-1e-10))
    res_best=meta.stage1_result; phi_best=meta.stage1_phi; best_idx=meta.stage1_idx;
    history=meta.stage1_history; meta.monotonicity_guard_triggered=true;
else
    meta.monotonicity_guard_triggered=false;
end
meta.refinement_runtime_s=toc(refine_timer);
meta.refinement_solver_stats=comm_solver_stats_delta(refine_solver_start,comm_solver_telemetry('snapshot'));
meta.total_exact_refine_evals=meta.mechanism_evals+meta.polish_evals;
meta.total_exact_algorithm_evals=meta.stage1_exact_evals+meta.total_exact_refine_evals;
meta.final_obj_snr=res_best.obj_snr; meta.final_p_comm=res_best.p_comm;
if meta.stage1_result.feasible && res_best.feasible && meta.stage1_result.obj_snr>0
    meta.refinement_gain_dB=10*log10(res_best.obj_snr/meta.stage1_result.obj_snr);
else
    meta.refinement_gain_dB=NaN;
end
end

function [best,best_phi,best_idx,used,improved]=evaluate_ranked_joint_candidates( ...
    ord,score,phi_list,idx_list,budget,current,Hd_u0,Hr_u,Hd_t0,Hr_t,G0, ...
    gain_u,gain_t,gain_ris,gamma,sigma_c,sigma_r2,P_total,omega,opts,skip)
if nargin<21; skip=0; end
best=current; best_phi=[]; best_idx=NaN; used=0; improved=false;
valid_seen=0;
for jj=1:numel(ord)
    q=ord(jj); if ~isfinite(score(q)); continue; end
    valid_seen=valid_seen+1; if valid_seen<=skip; continue; end
    if used>=max(budget-skip,0); break; end
    ti=idx_list(q);
    [Hu_t,Ht_t,G_t]=apply_ret_channels(ti,Hd_u0,Hd_t0,G0,gain_u,gain_t,gain_ris);
    rr=evaluate_dynamic_isac(phi_list{q},true,Hu_t,Hr_u,Ht_t,Hr_t,G_t, ...
        gamma,sigma_c,sigma_r2,P_total,omega);
    used=used+1;
    if exact_improvement(rr,best,opts)
        best=rr; best_phi=phi_list{q}; best_idx=ti; improved=true;
    end
end
end

function [phi,best_idx,res,history,meta]=sequential_multifidelity_trust_region_polish( ...
    phi,best_idx,res,history,Hd_u0,Hr_u,Hd_t0,Hr_t,G0,gain_u,gain_t,gain_ris, ...
    theta_range,gamma,sigma_c,sigma_r2,P_total,omega,opts)
%SEQUENTIAL_MULTIFIDELITY_TRUST_REGION_POLISH MA-RCG-MF v2 Stage-II.
%
% Diagnosis of the previous refinement showed that the dominant gap to the
% exact block-search reference came from restrictive candidate generation.
% V2 fixes that failure mode without turning the proposed method into an
% exhaustive search:
%   1) every RIS block is visited (no active-block hard pruning);
%   2) accepted block moves update the incumbent immediately, so improvements
%      from different blocks can accumulate within the same sweep;
%   3) the first sweep uses a full-circle *proxy-only* coarse probe, while
%      later sweeps use local candidates inside a per-block trust region;
%   4) only a small adaptive Top-K is exact-evaluated;
%   5) each block's trust radius expands after boundary improvements and
%      shrinks after unsuccessful/small local moves.
%
% Every accepted update is still checked by the exact communication-QoS
% solver and the exact structured sensing objective.

meta=struct('evals',0,'accepts',0,'screened',0,'round_accepts',[], ...
    'sweeps_completed',0,'block_visits',0,'budget_expansions',0, ...
    'ret_evals',0,'ris_evals',0,'sweep_gain_dB',[], ...
    'sweep_accepts',[],'trust_final_deg',[],'block_accept_count',[], ...
    'block_exact_evals',[],'accepted_offset_deg',[]);
if ~res.feasible; return; end

N=numel(phi);
B=max(1,min(N,round(get_field_or(opts,'v2_blocks',8))));
max_sweeps=max(1,round(get_field_or(opts,'v2_max_sweeps',3)));
trust_init=max(eps,real(get_field_or(opts,'v2_trust_init_deg',45)));
trust_min=max(eps,real(get_field_or(opts,'v2_trust_min_deg',2.8125)));
trust_max=max(trust_init,real(get_field_or(opts,'v2_trust_max_deg',135)));
trust_expand=max(1,real(get_field_or(opts,'v2_trust_expand',1.5)));
trust_shrink=min(max(real(get_field_or(opts,'v2_trust_shrink',0.5)),0.05),0.95);
q0=max(1,round(get_field_or(opts,'v2_exact_topk_initial',2)));
qmax=max(q0,round(get_field_or(opts,'v2_exact_topk_max',4)));
expand_budget=logical(get_field_or(opts,'v2_expand_budget_on_stall',true));
tie_gap_db=max(0,real(get_field_or(opts,'v2_proxy_tie_gap_dB',0.05)));
ret_half=max(0,round(get_field_or(opts,'v2_ret_half_window',1)));
ret_q0=max(1,round(get_field_or(opts,'v2_ret_exact_topk_initial',1)));
ret_qmax=max(ret_q0,round(get_field_or(opts,'v2_ret_exact_topk_max',2)));
max_noaccept=max(1,round(get_field_or(opts,'v2_max_noaccept_sweeps',2)));
sweep_stop_gain=max(0,real(get_field_or(opts,'v2_sweep_stop_gain_dB',0)));

trust=trust_init*ones(1,B);
block_id=ceil((1:N)*B/N);
block_accept_count=zeros(1,B);
block_exact_evals=zeros(1,B);
last_offset=nan(1,B);
noaccept_sweeps=0;

for sw=1:max_sweeps
    sweep_obj0=res.obj_snr;
    sweep_accepts=0;

    % ---- RET micro-search: cheap ranking, exact acceptance.  The accepted
    % RET point becomes the geometry used by all following RIS block visits.
    if ret_half>0
        ids=max(1,best_idx-ret_half):min(numel(theta_range),best_idx+ret_half);
        score=-Inf(size(ids));
        for jj=1:numel(ids)
            ti=ids(jj); if ti==best_idx; continue; end
            [Hu_t,Ht_t,G_t]=apply_ret_channels(ti,Hd_u0,Hd_t0,G0,gain_u,gain_t,gain_ris);
            score(jj)=mechanism_proxy(phi,Hu_t,Hr_u,Ht_t,Hr_t,G_t,res, ...
                gamma,sigma_c,sigma_r2,P_total,omega,opts);
            meta.screened=meta.screened+1;
        end
        [~,ord_ret]=sort(score,'descend');
        [cand,cidx,used,imp]=eval_v2_ret_ranked(ord_ret,score,ids,ret_q0,res,phi, ...
            Hd_u0,Hr_u,Hd_t0,Hr_t,G0,gain_u,gain_t,gain_ris,gamma,sigma_c, ...
            sigma_r2,P_total,omega,opts,0);
        meta.evals=meta.evals+used; meta.ret_evals=meta.ret_evals+used;
        if ~imp && expand_budget && ret_qmax>ret_q0
            [cand2,cidx2,used2,imp2]=eval_v2_ret_ranked(ord_ret,score,ids,ret_qmax,res,phi, ...
                Hd_u0,Hr_u,Hd_t0,Hr_t,G0,gain_u,gain_t,gain_ris,gamma,sigma_c, ...
                sigma_r2,P_total,omega,opts,ret_q0);
            meta.evals=meta.evals+used2; meta.ret_evals=meta.ret_evals+used2;
            if used2>0; meta.budget_expansions=meta.budget_expansions+1; end
            if imp2; cand=cand2; cidx=cidx2; imp=true; end
        end
        if imp && cidx~=best_idx && v2_exact_improvement(cand,res,opts)
            res=cand; best_idx=cidx; meta.accepts=meta.accepts+1;
            sweep_accepts=sweep_accepts+1;
            history=append_history(history,res.p_comm,res.obj_snr,theta_range(best_idx),NaN);
        end
    end

    [Hd_u,Hd_t,G]=apply_ret_channels(best_idx,Hd_u0,Hd_t0,G0,gain_u,gain_t,gain_ris);

    % Alternate forward/reverse sweeps to reduce a deterministic block-order
    % bias while still visiting every block exactly once per sweep.
    if mod(sw,2)==1; block_order=1:B; else; block_order=B:-1:1; end

    for oo=1:numel(block_order)
        bb=block_order(oo);
        elem=find(block_id==bb);
        if isempty(elem); continue; end
        meta.block_visits=meta.block_visits+1;

        offsets=v2_candidate_offsets(trust(bb),sw,opts);
        cand_phi=cell(1,numel(offsets));
        cand_score=-Inf(1,numel(offsets));
        for cc=1:numel(offsets)
            q=phi;
            q(elem)=q(elem)*exp(1i*offsets(cc)*pi/180);
            q=normalize_phase(q,N);
            cand_phi{cc}=q;
            cand_score(cc)=mechanism_proxy(q,Hd_u,Hr_u,Hd_t,Hr_t,G,res, ...
                gamma,sigma_c,sigma_r2,P_total,omega,opts);
            meta.screened=meta.screened+1;
        end
        [~,ord]=sort(cand_score,'descend');

        % If several top proxy scores are nearly tied, spend the larger exact
        % budget immediately because the diagnostic showed only moderate rank
        % correlation. Otherwise start with the small budget and expand only
        % after a stall.
        qfirst=q0;
        if v2_proxy_is_ambiguous(cand_score,ord,q0,qmax,tie_gap_db)
            qfirst=qmax; meta.budget_expansions=meta.budget_expansions+1;
        end
        [local_best,local_phi,best_off,used,imp]=eval_v2_block_ranked( ...
            ord,cand_score,cand_phi,offsets,qfirst,res,Hd_u,Hr_u,Hd_t,Hr_t,G, ...
            gamma,sigma_c,sigma_r2,P_total,omega,opts,0);
        meta.evals=meta.evals+used; meta.ris_evals=meta.ris_evals+used;
        block_exact_evals(bb)=block_exact_evals(bb)+used;

        if ~imp && expand_budget && qfirst<qmax
            [best2,phi2,off2,used2,imp2]=eval_v2_block_ranked( ...
                ord,cand_score,cand_phi,offsets,qmax,res,Hd_u,Hr_u,Hd_t,Hr_t,G, ...
                gamma,sigma_c,sigma_r2,P_total,omega,opts,qfirst);
            meta.evals=meta.evals+used2; meta.ris_evals=meta.ris_evals+used2;
            block_exact_evals(bb)=block_exact_evals(bb)+used2;
            if used2>0; meta.budget_expansions=meta.budget_expansions+1; end
            if imp2; local_best=best2; local_phi=phi2; best_off=off2; imp=true; end
        end

        if imp && v2_exact_improvement(local_best,res,opts)
            % Sequential update: unlike R3, the next block sees this newly
            % accepted incumbent, allowing cooperative multi-block gains.
            old_obj=res.obj_snr;
            res=local_best; phi=local_phi;
            meta.accepts=meta.accepts+1; sweep_accepts=sweep_accepts+1;
            block_accept_count(bb)=block_accept_count(bb)+1;
            last_offset(bb)=best_off;
            history=append_history(history,res.p_comm,res.obj_snr,theta_range(best_idx),NaN);

            abs_off=abs(wrap_offset_deg(best_off));
            old_r=trust(bb);
            if abs_off>old_r*(1+1e-9) || abs_off>=0.75*old_r
                trust(bb)=min(trust_max,max(old_r*trust_expand,abs_off));
            elseif abs_off<=0.35*old_r
                % A successful small step indicates that the useful local
                % scale is narrowing; contract gently rather than abruptly.
                trust(bb)=max(trust_min,0.75*old_r);
            end
            if ~(res.obj_snr>old_obj); error('V2 monotonicity violation inside block update.'); end
        else
            trust(bb)=max(trust_min,trust(bb)*trust_shrink);
        end
    end

    meta.sweeps_completed=sw;
    meta.sweep_accepts(sw)=sweep_accepts; %#ok<AGROW>
    if sweep_obj0>0 && res.obj_snr>0
        meta.sweep_gain_dB(sw)=10*log10(res.obj_snr/sweep_obj0); %#ok<AGROW>
    else
        meta.sweep_gain_dB(sw)=NaN; %#ok<AGROW>
    end
    meta.round_accepts(sw)=sweep_accepts>0; %#ok<AGROW>

    if sweep_accepts==0
        noaccept_sweeps=noaccept_sweeps+1;
    else
        noaccept_sweeps=0;
    end

    % A no-accept sweep is not immediately treated as convergence: trust
    % regions have already contracted and are allowed to retry. Stop only
    % after repeated no-accept sweeps or after all regions reach the floor.
    if noaccept_sweeps>=max_noaccept || ...
            (sweep_accepts==0 && all(trust<=trust_min*(1+1e-12)))
        break;
    end
    % If a late sweep produces only a numerically tiny total gain and every
    % trust region is already local, further exact solves are not justified.
    if sw>=2 && sweep_accepts>0 && isfinite(meta.sweep_gain_dB(sw)) && ...
            meta.sweep_gain_dB(sw)<sweep_stop_gain && ...
            all(trust<=max(2*trust_min,trust_init/2))
        break;
    end
end

meta.trust_final_deg=trust;
meta.block_accept_count=block_accept_count;
meta.block_exact_evals=block_exact_evals;
meta.accepted_offset_deg=last_offset;
end

function offsets=v2_candidate_offsets(radius_deg,sweep,opts)
% Local trust-region offsets plus a proxy-only global coarse grid on the
% earliest sweep(s).  All angles are represented on the shortest [-180,180]
% arc and duplicate offsets are removed.
frac=real(get_field_or(opts,'v2_local_fractions',[1 0.5 0.25]));
frac=unique(abs(frac(isfinite(frac)&frac>0)),'stable');
if isempty(frac); frac=[1 0.5 0.25]; end
local=[];
for ff=1:numel(frac)
    d=radius_deg*frac(ff); local=[local -d d]; %#ok<AGROW>
end
offsets=local;
nglobal=max(0,round(get_field_or(opts,'v2_global_probe_levels',0)));
gsweeps=max(0,round(get_field_or(opts,'v2_global_probe_sweeps',1)));
if sweep<=gsweeps && nglobal>=4
    g=360*(1:nglobal-1)/nglobal;
    g=arrayfun(@wrap_offset_deg,g);
    offsets=[offsets g]; %#ok<AGROW>
end
offsets=arrayfun(@wrap_offset_deg,offsets);
offsets=offsets(abs(offsets)>1e-10);
% Stable tolerance-aware uniqueness.
key=round(offsets*1e8)/1e8;
[~,ia]=unique(key,'stable'); offsets=offsets(sort(ia));
end

function tf=v2_proxy_is_ambiguous(score,ord,q0,qmax,tie_gap_db)
tf=false;
if qmax<=q0 || isempty(ord); return; end
ids=ord(isfinite(score(ord)) & score(ord)>0);
if numel(ids)<min(qmax,2); return; end
k=min(qmax,numel(ids));
s1=score(ids(1)); sk=score(ids(k));
if s1<=0 || sk<=0; return; end
gap_db=10*log10(s1/sk);
tf=isfinite(gap_db) && gap_db<=tie_gap_db;
end

function [best,bphi,best_off,used,imp]=eval_v2_block_ranked( ...
    ord,score,phi_list,offsets,budget,current,Hd_u,Hr_u,Hd_t,Hr_t,G, ...
    gamma,sigma_c,sigma_r2,P_total,omega,opts,skip)
if nargin<18; skip=0; end
best=current; bphi=[]; best_off=NaN; used=0; imp=false; seen=0;
for jj=1:numel(ord)
    q=ord(jj); if ~isfinite(score(q)); continue; end
    seen=seen+1; if seen<=skip; continue; end
    if used>=max(budget-skip,0); break; end
    x=evaluate_dynamic_isac(phi_list{q},true,Hd_u,Hr_u,Hd_t,Hr_t,G, ...
        gamma,sigma_c,sigma_r2,P_total,omega);
    used=used+1;
    if v2_exact_improvement(x,best,opts)
        best=x; bphi=phi_list{q}; best_off=offsets(q); imp=true;
    end
end
end

function [cand,cidx,used,imp]=eval_v2_ret_ranked(ord,score,ids,budget,current,phi, ...
    Hd_u0,Hr_u,Hd_t0,Hr_t,G0,gain_u,gain_t,gain_ris,gamma,sigma_c,sigma_r2, ...
    P_total,omega,opts,skip)
if nargin<21; skip=0; end
cand=current; cidx=NaN; used=0; imp=false; seen=0;
for jj=1:numel(ord)
    q=ord(jj); if ~isfinite(score(q)); continue; end
    seen=seen+1; if seen<=skip; continue; end
    if used>=max(budget-skip,0); break; end
    ti=ids(q);
    [Hu_t,Ht_t,G_t]=apply_ret_channels(ti,Hd_u0,Hd_t0,G0,gain_u,gain_t,gain_ris);
    x=evaluate_dynamic_isac(phi,true,Hu_t,Hr_u,Ht_t,Hr_t,G_t, ...
        gamma,sigma_c,sigma_r2,P_total,omega);
    used=used+1;
    if v2_exact_improvement(x,cand,opts); cand=x; cidx=ti; imp=true; end
end
end

function tf=v2_exact_improvement(candidate,current,opts)
tf=false;
if ~isstruct(candidate)||~isfield(candidate,'feasible')||~candidate.feasible; return; end
if ~isstruct(current)||~isfield(current,'feasible')||~current.feasible; tf=true; return; end
if ~isfinite(candidate.obj_snr)||~isfinite(current.obj_snr)||current.obj_snr<=0; return; end
min_gain_db=max(0,real(get_field_or(opts,'v2_min_gain_dB',0)));
ratio=max(1+get_field_or(opts,'refine_accept_tol',1e-10),10^(min_gain_db/10));
tf=candidate.obj_snr>current.obj_snr*ratio;
end

function d=wrap_offset_deg(d)
d=mod(d+180,360)-180;
% Keep +180 rather than -180 only for display/consistency; both are equal.
if abs(d+180)<1e-10; d=180; end
end

function [phi,best_idx,res,history,meta]=adaptive_multifidelity_polish( ...
    phi,best_idx,res,history,Hd_u0,Hr_u,Hd_t0,Hr_t,G0,gain_u,gain_t,gain_ris, ...
    theta_range,gamma,sigma_c,sigma_r2,P_total,omega,opts)
meta=struct('evals',0,'accepts',0,'screened',0,'round_accepts',[]);
if ~res.feasible; return; end
Bsch=get_field_or(opts,'polish_blocks_schedule',4); Bsch=max(1,round(Bsch(:).'));
stepSch=get_field_or(opts,'polish_phase_step_schedule',{11.25});
if ~iscell(stepSch); stepSch={stepSch}; end
Asch=get_field_or(opts,'polish_active_blocks_schedule',Bsch); Asch=max(1,round(Asch(:).'));
Q0sch=get_field_or(opts,'polish_exact_topk_initial_schedule',2); Q0sch=max(1,round(Q0sch(:).'));
Qmaxsch=get_field_or(opts,'polish_exact_topk_max_schedule',Q0sch); Qmaxsch=max(1,round(Qmaxsch(:).'));
R=max([numel(Bsch),numel(stepSch),numel(Asch),numel(Q0sch),numel(Qmaxsch)]);
ret_half=max(0,round(get_field_or(opts,'polish_ret_half_window',1)));
ret_q0=max(1,round(get_field_or(opts,'polish_ret_exact_topk_initial',1)));
ret_qmax=max(ret_q0,round(get_field_or(opts,'polish_ret_exact_topk_max',2)));
try_fine=logical(get_field_or(opts,'polish_try_fine_on_stall',true));
N=numel(phi);
stall_rounds=0;

for rrnd=1:R
    round_improved=false;
    B=min(N,get_sched(Bsch,rrnd)); Aactive=min(B,get_sched(Asch,rrnd));
    steps=abs(get_sched_cell(stepSch,rrnd)); steps=steps(isfinite(steps)&steps>0);
    if isempty(steps); steps=11.25; end
    q0=get_sched(Q0sch,rrnd); qmax=max(q0,get_sched(Qmaxsch,rrnd));

    % RET micro-search with adaptive exact budget.
    if ret_half>0
        ids=max(1,best_idx-ret_half):min(numel(theta_range),best_idx+ret_half);
        score=-Inf(size(ids));
        for jj=1:numel(ids)
            ti=ids(jj); if ti==best_idx; continue; end
            [Hu_t,Ht_t,G_t]=apply_ret_channels(ti,Hd_u0,Hd_t0,G0,gain_u,gain_t,gain_ris);
            score(jj)=mechanism_proxy(phi,Hu_t,Hr_u,Ht_t,Hr_t,G_t,res, ...
                gamma,sigma_c,sigma_r2,P_total,omega,opts);
            meta.screened=meta.screened+1;
        end
        [~,ord]=sort(score,'descend');
        [cand,cidx,used,imp]=eval_ret_ranked(ord,score,ids,ret_q0,res,phi, ...
            Hd_u0,Hr_u,Hd_t0,Hr_t,G0,gain_u,gain_t,gain_ris,gamma,sigma_c, ...
            sigma_r2,P_total,omega,opts,0);
        meta.evals=meta.evals+used;
        if ~imp && ret_qmax>ret_q0
            [cand2,cidx2,used2,imp2]=eval_ret_ranked(ord,score,ids,ret_qmax,res,phi, ...
                Hd_u0,Hr_u,Hd_t0,Hr_t,G0,gain_u,gain_t,gain_ris,gamma,sigma_c, ...
                sigma_r2,P_total,omega,opts,ret_q0);
            meta.evals=meta.evals+used2;
            if imp2; cand=cand2; cidx=cidx2; imp=true; end
        end
        if imp && cidx~=best_idx && exact_improvement(cand,res,opts)
            res=cand; best_idx=cidx; meta.accepts=meta.accepts+1; round_improved=true;
            history=append_history(history,res.p_comm,res.obj_snr,theta_range(best_idx),NaN);
        end
    end

    [Hd_u,Hd_t,G]=apply_ret_channels(best_idx,Hd_u0,Hd_t0,G0,gain_u,gain_t,gain_ris);
    block_id=ceil((1:N)*B/N);
    cand_phi=cell(0,1); cand_score=[]; cand_block=[];
    for bb=1:B
        elem=find(block_id==bb);
        for ss=1:numel(steps)
            for sg=[-1 1]
                q=phi; q(elem)=q(elem)*exp(1i*sg*steps(ss)*pi/180);
                q=normalize_phase(q,N);
                sc=mechanism_proxy(q,Hd_u,Hr_u,Hd_t,Hr_t,G,res, ...
                    gamma,sigma_c,sigma_r2,P_total,omega,opts);
                cand_phi{end+1}=q; cand_score(end+1)=sc; cand_block(end+1)=bb; %#ok<AGROW>
                meta.screened=meta.screened+1;
            end
        end
    end

    % Active-block selection: rank blocks by their best proxy candidate.
    block_best=-Inf(1,B);
    for bb=1:B
        z=cand_score(cand_block==bb); if ~isempty(z); block_best(bb)=max(z); end
    end
    [~,bOrd]=sort(block_best,'descend'); active=bOrd(1:min(Aactive,sum(isfinite(block_best))));
    active_mask=ismember(cand_block,active)&isfinite(cand_score);
    ord=diverse_candidate_order(cand_score,cand_block,active_mask,active);

    [local_best,local_phi,used,imp]=eval_ris_ranked(ord,cand_score,cand_phi,q0,res, ...
        Hd_u,Hr_u,Hd_t,Hr_t,G,gamma,sigma_c,sigma_r2,P_total,omega,opts,0);
    meta.evals=meta.evals+used;
    if ~imp && qmax>q0
        [best2,phi2,used2,imp2]=eval_ris_ranked(ord,cand_score,cand_phi,qmax,res, ...
            Hd_u,Hr_u,Hd_t,Hr_t,G,gamma,sigma_c,sigma_r2,P_total,omega,opts,q0);
        meta.evals=meta.evals+used2;
        if imp2; local_best=best2; local_phi=phi2; imp=true; end
    end
    if imp && exact_improvement(local_best,res,opts)
        res=local_best; phi=local_phi; meta.accepts=meta.accepts+1; round_improved=true;
        history=append_history(history,res.p_comm,res.obj_snr,theta_range(best_idx),NaN);
    end
    meta.round_accepts(rrnd)=round_improved; %#ok<AGROW>
    if round_improved; stall_rounds=0; else; stall_rounds=stall_rounds+1; end
    if ~round_improved && (~try_fine || stall_rounds>=2 || rrnd>=R); break; end
end
end

function [cand,cidx,used,imp]=eval_ret_ranked(ord,score,ids,budget,current,phi, ...
    Hd_u0,Hr_u,Hd_t0,Hr_t,G0,gain_u,gain_t,gain_ris,gamma,sigma_c,sigma_r2, ...
    P_total,omega,opts,skip)
cand=current; cidx=NaN; used=0; imp=false; seen=0;
for jj=1:numel(ord)
    q=ord(jj); if ~isfinite(score(q)); continue; end
    seen=seen+1; if seen<=skip; continue; end
    if used>=max(budget-skip,0); break; end
    ti=ids(q); [Hu_t,Ht_t,G_t]=apply_ret_channels(ti,Hd_u0,Hd_t0,G0,gain_u,gain_t,gain_ris);
    x=evaluate_dynamic_isac(phi,true,Hu_t,Hr_u,Ht_t,Hr_t,G_t,gamma,sigma_c,sigma_r2,P_total,omega);
    used=used+1; if exact_improvement(x,cand,opts); cand=x; cidx=ti; imp=true; end
end
end

function [best,bphi,used,imp]=eval_ris_ranked(ord,score,phi_list,budget,current, ...
    Hd_u,Hr_u,Hd_t,Hr_t,G,gamma,sigma_c,sigma_r2,P_total,omega,opts,skip)
best=current; bphi=[]; used=0; imp=false; seen=0;
for jj=1:numel(ord)
    q=ord(jj); if ~isfinite(score(q)); continue; end
    seen=seen+1; if seen<=skip; continue; end
    if used>=max(budget-skip,0); break; end
    x=evaluate_dynamic_isac(phi_list{q},true,Hd_u,Hr_u,Hd_t,Hr_t,G, ...
        gamma,sigma_c,sigma_r2,P_total,omega);
    used=used+1; if exact_improvement(x,best,opts); best=x; bphi=phi_list{q}; imp=true; end
end
end

function ord=diverse_candidate_order(score,block,mask,active)
% First take the strongest candidate from each active block, then append the
% remaining active candidates globally. This prevents the proxy top-K from
% being monopolized by one block when its local surrogate is slightly biased.
ord=[]; chosen=false(size(score));
for bb=active(:).'
    ids=find(mask & block==bb); if isempty(ids); continue; end
    [~,j]=max(score(ids)); q=ids(j); ord(end+1)=q; chosen(q)=true; %#ok<AGROW>
end
if ~isempty(ord)
    [~,z]=sort(score(ord),'descend'); ord=ord(z);
end
rest=find(mask & ~chosen);
[~,z]=sort(score(rest),'descend'); ord=[ord rest(z)];
end

function v=get_sched(x,i)
x=x(:).'; v=x(min(i,numel(x)));
end

function v=get_sched_cell(c,i)
v=c{min(i,numel(c))};
end

function score=mechanism_proxy(phi,Hd_u,Hr_u,Hd_t,Hr_t,G,current, ...
    gamma,sigma_c,sigma_r2,P_total,omega,opts) %#ok<INUSD>
%MECHANISM_PROXY Cheap power--channel coupled ranking metric.
%
% For a candidate passive configuration, keep the *directions* of the current
% communication beamformers and compute the minimum common power-scaling
% factor a_req that makes those fixed directions satisfy every SINR target:
%
%   a_req |h_k^H w_k|^2 >= gamma_k
%       [a_req sum_{j~=k}|h_k^H w_j|^2 + sigma_c^2].
%
% Hence
%   a_req = max_k gamma_k sigma_c^2 /
%                    (S_k - gamma_k I_k),
% whenever every denominator is positive.  This gives a feasible (although
% generally not minimum-power) communication design for the candidate, so
% Pcomm_hat=a_req||Wc||_F^2 is an inexpensive upper estimate of its reoptimized
% minimum communication power.  The structured sensing decomposition then
% gives the screening objective
%
%   J_proxy = a_req J_comm,fixed + (P_T-Pcomm_hat) G_eff.
%
% The proxy is used ONLY for ranking local candidates.  Final acceptance is
% always based on a validated exact communication-QoS solve (duality FPI
% with automatic CVX fallback) and the exact structured radar objective.
Hu=build_effective_user_channel(phi,true,Hd_u,Hr_u,G);
C=build_radar_matrix(true,phi,Hd_t,Hr_t,G,sigma_r2,omega);
Geff=effective_sensing_gain(Hu,C);
a_req=uniform_comm_power_scale(Hu,current.Wc,gamma,sigma_c);
if ~isfinite(a_req); score=-Inf; return; end
Pbase=real(norm(current.Wc,'fro')^2);
Pcomm_hat=a_req*Pbase;
if ~isfinite(Pcomm_hat) || Pcomm_hat>P_total*(1+1e-9)
    score=-Inf; return;
end
Jcomm_fixed=real(trace(current.Wc'*C*current.Wc));
Jcomm_hat=max(a_req,0)*max(Jcomm_fixed,0);
Pr_hat=max(P_total-Pcomm_hat,0);
score=max(real(Jcomm_hat+Pr_hat*Geff),0);
if ~isfinite(score); score=-Inf; end
end

function a_req=uniform_comm_power_scale(Hu,Wc,gamma,sigma_c)
% Minimum common power scale for fixed communication beam directions.
K=size(Hu,2); noise=sigma_c^2; a_req=0;
if isscalar(gamma); gamma_vec=repmat(gamma,K,1); else; gamma_vec=gamma(:); end
if numel(gamma_vec)~=K; a_req=Inf; return; end
for k=1:K
    hk=Hu(:,k);
    sig=abs(hk'*Wc(:,k))^2;
    interf=0;
    for j=1:size(Wc,2)
        if j~=k; interf=interf+abs(hk'*Wc(:,j))^2; end
    end
    margin=sig-gamma_vec(k)*interf;
    if ~isfinite(margin) || margin<=1e-18
        a_req=Inf; return;
    end
    ak=gamma_vec(k)*noise/margin;
    if ~isfinite(ak) || ak<0; a_req=Inf; return; end
    a_req=max(a_req,ak);
end
end

function tf=exact_improvement(candidate,current,opts)
tf=false;
if ~isstruct(candidate) || ~isfield(candidate,'feasible') || ~candidate.feasible
    return;
end
if ~isstruct(current) || ~isfield(current,'feasible') || ~current.feasible
    tf=true; return;
end
if ~isfinite(candidate.obj_snr) || ~isfinite(current.obj_snr) || current.obj_snr<=0
    return;
end
min_gain_db=get_field_or(opts,'polish_min_gain_dB',0);
ratio=max(1+get_field_or(opts,'refine_accept_tol',1e-10),10^(min_gain_db/10));
tf=candidate.obj_snr>current.obj_snr*ratio;
end

function phi=phase_geodesic_blend(phi_a,phi_b,alpha)
% Move alpha of the shortest element-wise phase arc from phi_a to phi_b.
N=numel(phi_a); a=normalize_phase(phi_a,N); b=normalize_phase(phi_b,N);
alpha=min(max(real(alpha),0),1);
delta=angle(b.*conj(a));
phi=a.*exp(1i*alpha*delta);
phi=normalize_phase(phi,N);
end

function [phi_best,res_best]=select_phase_radar(phi_warm,num_trials, ...
    Hd_u,Hr_u,Hd_t,Hr_t,G,gamma,sigma_c,sigma_r2,P_total,omega)
N=size(G,1); T=size(Hd_t,2);
C=cell(0,1);
C{end+1}=normalize_phase(phi_warm,N);
C{end+1}=target_focus_phase(Hd_t,Hr_t,G,omega,[]);
C{end+1}=communication_focus_phase(Hd_u,Hr_u,G);
C{end+1}=exp(1i*2*pi*(0:N-1)'/N);
center_t=ceil(T/2);
C{end+1}=target_focus_phase(Hd_t,Hr_t,G,omega,center_t);
C{end+1}=exp(1i*2*pi*rand(N,1));
while numel(C)<num_trials; C{end+1}=exp(1i*2*pi*rand(N,1)); end
C=C(1:min(numel(C),max(num_trials,4)));

phi_best=C{1}; res_best=empty_result(); res_best.obj_snr=-Inf;
for q=1:numel(C)
    rr=evaluate_dynamic_isac(C{q},true,Hd_u,Hr_u,Hd_t,Hr_t,G, ...
        gamma,sigma_c,sigma_r2,P_total,omega);
    if rr.feasible && rr.obj_snr>res_best.obj_snr
        res_best=rr; phi_best=C{q};
    end
end
end

%% ======================== Evaluation =====================================
function comm=evaluate_comm(phi,use_ris,Hd_u,Hr_u,G,gamma,sigma_c,P_total)
comm=empty_comm();
if use_ris
    Hu=build_effective_user_channel(phi,true,Hd_u,Hr_u,G);
else
    Hu=Hd_u;
end
[Wc,p_comm,feas]=solve_comm_min_power(Hu,gamma,sigma_c);
if ~feas || ~isfinite(p_comm) || p_comm>P_total*(1+1e-6); return; end
comm.feasible=true; comm.p_comm=p_comm; comm.Wc=Wc; comm.Hu=Hu;
end

function res=evaluate_dynamic_isac(phi,use_ris,Hd_u,Hr_u,Hd_t,Hr_t,G, ...
    gamma,sigma_c,sigma_r2,P_total,omega)
comm=evaluate_comm(phi,use_ris,Hd_u,Hr_u,G,gamma,sigma_c,P_total);
res=finalize_isac_from_comm(comm,phi,use_ris,Hd_t,Hr_t,G, ...
    sigma_r2,P_total,omega);
end

function res=finalize_isac_from_comm(comm,phi,use_ris,Hd_t,Hr_t,G, ...
    sigma_r2,P_total,omega)
res=empty_result();
if ~comm.feasible; return; end
C=build_radar_matrix(use_ris,phi,Hd_t,Hr_t,G,sigma_r2,omega);
[Wr,p_radar]=allocate_radar_power(comm.Hu,C,P_total,comm.p_comm);
W=[comm.Wc,Wr];
res.feasible=true; res.p_comm=comm.p_comm; res.p_radar=p_radar;
res.obj_snr=real(trace(W'*C*W)); res.Wc=comm.Wc; res.Wr=Wr;
res.Hu=comm.Hu; res.C=C;
end

%% ======================== Initializers / Wrappers =========================
function phi=communication_focus_phase(Hd_u,Hr_u,G)
N=size(G,1); score=zeros(N,1); K=size(Hd_u,2);
for k=1:K
    gk=conj(G)*conj(Hd_u(:,k));
    score=score+2*conj(Hr_u(:,k).*gk)/K;
end
phi=phase_from_score(score,N);
end

function phi=target_focus_phase(Hd_t,Hr_t,G,omega,target_idx)
N=size(G,1); score=zeros(N,1); T=size(Hd_t,2);
if nargin<5 || isempty(target_idx); ids=1:T; else; ids=target_idx; end
for t=ids
    gt=conj(G)*conj(Hd_t(:,t));
    f2=2*conj(Hr_t(:,t).*gt);
    if numel(ids)==1 && ~isempty(target_idx); wt=1; else; wt=omega(t); end
    score=score+wt*f2;
end
phi=phase_from_score(score,N);
end

function phi=phase_from_score(score,N)
if isempty(score) || norm(score)<1e-14
    phi=ones(N,1);
else
    phi=exp(-1i*angle(score(:)));
end
phi=normalize_phase(phi,N);
end

function phi=call_opt_phi_minpc_warm(N,K,M,Hd_u,Hr_u,G,W,A,rho,phi0)
try
    phi=opt_phi_min_pc(N,K,M,Hd_u,Hr_u,G,W,A,rho,phi0);
catch ME
    if contains(lower(ME.message),'too many input')
        phi=opt_phi_min_pc(N,K,M,Hd_u,Hr_u,G,W,A,rho);
    else
        rethrow(ME);
    end
end
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

%% ======================== RET / Numeric Helpers ===========================
function [Hd_u,Hd_t,Gt]=apply_ret_channels(ti,Hd_u0,Hd_t0,G0,gain_u,gain_t,gain_ris)
Hd_u=Hd_u0; Hd_t=Hd_t0;
for k=1:size(Hd_u0,2); Hd_u(:,k)=sqrt(gain_u(ti,k))*Hd_u0(:,k); end
for t=1:size(Hd_t0,2); Hd_t(:,t)=sqrt(gain_t(ti,t))*Hd_t0(:,t); end
Gt=sqrt(gain_ris(ti))*G0;
end

function [Hd_u,Gt]=apply_ret_comm_channels(ti,Hd_u0,G0,gain_u,gain_ris)
Hd_u=Hd_u0;
for k=1:size(Hd_u0,2); Hd_u(:,k)=sqrt(gain_u(ti,k))*Hd_u0(:,k); end
Gt=sqrt(gain_ris(ti))*G0;
end

function A=make_A(Hu,Wc)
K=size(Hu,2); A=zeros(K,K);
for k=1:K
    for j=1:K; A(k,j)=Hu(:,k)'*Wc(:,j); end
end
end

function r=auxiliary_residual(A,Hu,Wc)
Aphys=make_A(Hu,Wc); d=A(:)-Aphys(:); r=max(abs(d));
end

function phi=normalize_phase(phi,N)
if isempty(phi) || numel(phi)~=N
    phi=ones(N,1);
else
    phi=phi(:);
end
bad=abs(phi)<1e-12; phi(bad)=1; phi=phi./abs(phi);
end

function h=empty_history()
h=struct('Pcomm',[],'ObjSNR',[],'ThetaDeg',[],'AuxResidual',[]);
end

function h=append_history(h,pcomm,obj,theta,aux)
h.Pcomm(end+1)=pcomm;
h.ObjSNR(end+1)=obj;
h.ThetaDeg(end+1)=theta;
h.AuxResidual(end+1)=aux;
end

function h=single_point_history(res,theta)
h=empty_history();
if isfield(res,'feasible') && res.feasible
    h=append_history(h,res.p_comm,res.obj_snr,theta,NaN);
end
end

function comm=empty_comm()
comm=struct('feasible',false,'p_comm',Inf,'Wc',[],'Hu',[]);
end

function res=empty_result()
res=struct('feasible',false,'p_comm',NaN,'p_radar',NaN,'obj_snr',NaN, ...
    'Wc',[],'Wr',[],'Hu',[],'C',[]);
end
