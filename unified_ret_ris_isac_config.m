function cfg = unified_ret_ris_isac_config(run_mode)
%UNIFIED_RET_RIS_ISAC_CONFIG Final clean parameters for Experiments 1--12 (paper experiments).
% All experiments use the same antennas, roadside geometry, user/target
% model, path-loss assumptions, noise powers, RET model and random seeds.
%
% Exp. 1: scan communication SINR and minimize communication power.
% Exp. 2: fix communication SINR, scan total power and evaluate sensing.
% Exp. 3: fix total power, scan communication SINR and evaluate coupling.
% Exp. 4: RET/RIS ablation at a representative operating point.
% Exp. 5: RIS-size/user-load scalability; optional supplementary N/K timing.
% Exp. 6: multi-timescale mobility and stale RET-RIS configuration loss.
% Exp. 7: algorithmic baselines and true-objective convergence.
% Exp. 8: RIS roadside deployment sensitivity.
% Exp. 9: imperfect-CSI sensitivity on the true channel.
% Exp.10: post-quantization RIS hardware loss.
% Exp.11: RET angular-resolution/runtime tradeoff.
% Exp.12: complexity and paired statistical significance.

if nargin < 1 || isempty(run_mode)
    run_mode = 'fast';
end
run_mode = lower(char(run_mode));

cfg = struct();
cfg.RUN_MODE = run_mode;
cfg.paper_release_version = '2026-08-08-ma-rcg-mf-v2-sequential-mf-atr';
cfg.exp6_pairing_revision = '2026-08-08-exp6-strict-paired-crn-v1';

switch run_mode
    case 'fast'
        cfg.ITER = 3;
        cfg.gamma_dB_range_exp1 = 0:4:16;
        cfg.gamma_dB_range_exp3 = 0:4:16;
        cfg.max_ao_iter = 6;
        cfg.theta_tilt_range = 0:1:35; % same 1-deg hardware grid as PAPER mode
        cfg.tilt_local_half_window = 1;
        cfg.num_phase_trials = 6;
        cfg.num_workers = 3;
        cfg.ITER_EXP4 = 3;
        cfg.ITER_EXP5 = 5;
        cfg.N_range_exp5 = [32 64 96];
        cfg.K_range_exp5 = [2 6 8 10 12];
        % Experiment-5 numerical-quality and timing controls.
        cfg.exp5_execution_profile = 'fast';
        cfg.exp5_quality_restarts = 1;
        cfg.exp5_use_continuation = true;
        cfg.exp5_run_timing_benchmark = false;
        cfg.exp5_timing_mc = 1;
        cfg.exp5_timing_repeats = 1;
        cfg.exp5_timing_warmup_runs = 1;
        cfg.ITER_EXP6 = 6;
        cfg.exp6_bootstrap_reps = 250;
        cfg.speed_kmh_exp6 = [30 90 120];
        cfg.ris_update_ms_exp6 = [0 5 20 100];
        cfg.ITER_EXP7 = 4;
        cfg.exp7_random_trials = 12;
        cfg.exp7_bcs_blocks = 4;
        cfg.exp7_bcs_phase_levels = 8;
        cfg.exp7_bcs_sweeps = 2;
        % Fast-good optimizer: RCG + adaptive multi-fidelity exact evaluation.
        cfg.proposed_mechanism_fallback = true;
        cfg.proposed_stage1_geodesic_fractions = [1 0.5];
        cfg.proposed_stage1_exact_topk_initial = 1;
        cfg.proposed_stage1_exact_topk_max = 3;
        cfg.proposed_stage1_expand_on_stall = true;
        cfg.proposed_polish_blocks_schedule = [4 8];
        cfg.proposed_polish_phase_step_schedule = {[22.5 11.25],[11.25 5.625]};
        cfg.proposed_polish_active_blocks_schedule = [3 4];
        cfg.proposed_polish_exact_topk_initial_schedule = [2 2];
        cfg.proposed_polish_exact_topk_max_schedule = [3 4];
        cfg.proposed_polish_ret_half_window = 1;
        cfg.proposed_polish_ret_exact_topk_initial = 1;
        cfg.proposed_polish_ret_exact_topk_max = 2;
        cfg.proposed_polish_try_fine_on_stall = true;
        cfg.proposed_polish_min_gain_dB = 1e-4;
        % MA-RCG-MF v2 Stage-II: sequential multi-fidelity block refinement
        % with per-block adaptive trust regions. The legacy polish fields
        % above are retained for fair fixed-RIS baseline refinement only.
        cfg.proposed_stage2_variant = 'sequential_mf_atr';
        cfg.v2_blocks = 4;
        cfg.v2_max_sweeps = 2;
        cfg.v2_global_probe_levels = 8;
        cfg.v2_global_probe_sweeps = 1;
        cfg.v2_trust_init_deg = 45;
        cfg.v2_trust_min_deg = 5.625;
        cfg.v2_trust_max_deg = 135;
        cfg.v2_trust_expand = 1.5;
        cfg.v2_trust_shrink = 0.5;
        cfg.v2_local_fractions = [1 0.5 0.25];
        cfg.v2_exact_topk_initial = 1;
        cfg.v2_exact_topk_max = 3;
        cfg.v2_expand_budget_on_stall = true;
        cfg.v2_proxy_tie_gap_dB = 0.05;
        cfg.v2_ret_half_window = 1;
        cfg.v2_ret_exact_topk_initial = 1;
        cfg.v2_ret_exact_topk_max = 2;
        cfg.v2_min_gain_dB = 1e-3;
        cfg.v2_sweep_stop_gain_dB = 2e-3;
        cfg.v2_max_noaccept_sweeps = 2;
        % Engineering validation experiments.
        cfg.ITER_EXP8 = 3;
        cfg.ris_x_range_exp8 = [50 90 120];
        cfg.ITER_EXP9 = 5;
        cfg.exp9_bootstrap_reps = 300;
        cfg.csi_nmse_dB_range_exp9 = [-Inf -20 -10 -5];
        cfg.csi_error_range_exp9 = 10.^(cfg.csi_nmse_dB_range_exp9/20);
        cfg.ITER_EXP10 = 3;
        cfg.ris_phase_bits_exp10 = [1 2 3 Inf];
        cfg.ITER_EXP11 = 3;
        cfg.ret_resolution_deg_exp11 = [5 2 1 0.5];
        cfg.ITER_EXP12 = 3;
    case 'paper'
        cfg.ITER = 50;
        cfg.gamma_dB_range_exp1 = 0:2:16;
        cfg.gamma_dB_range_exp3 = 0:2:16;
        cfg.max_ao_iter = 10;
        cfg.theta_tilt_range = 0:1:35;
        cfg.tilt_local_half_window = 2;
        cfg.num_phase_trials = 10;
        % CVX jobs are memory intensive. Increase only after checking RAM.
        cfg.num_workers = 3;
        cfg.ITER_EXP4 = 50;
        cfg.ITER_EXP5 = 60;
        cfg.N_range_exp5 = [16 32 64 96 128];
        cfg.K_range_exp5 = [2 4 6 8 9 10 11 12];
        % Two exact-objective restarts reduce dimension-dependent local-optimum
        % bias. The expensive separate serial timing benchmark is disabled by
        % default because it repeats the complete N/K scans after performance
        % simulation. Exp.7/12 provide the manuscript runtime comparison.
        % Set RET_RIS_EXP5_FORCE_TIMING=1 only for an optional supplement.
        cfg.exp5_execution_profile = 'balanced';
        cfg.exp5_quality_restarts = 2;
        cfg.exp5_use_continuation = true;
        cfg.exp5_run_timing_benchmark = false;
        % Optional supplementary timing controls. These are ignored unless
        % RET_RIS_EXP5_FORCE_TIMING=1 is explicitly set. One serial sweep is
        % the recommended lightweight setting; full legacy timing remains
        % available only when an explicit supplementary timing profile is selected.
        cfg.exp5_timing_mc = 1;
        cfg.exp5_timing_repeats = 1;
        cfg.exp5_timing_warmup_runs = 1;
        cfg.ITER_EXP6 = 200;
        cfg.exp6_bootstrap_reps = 2000;
        cfg.speed_kmh_exp6 = [30 60 90 120];
        cfg.ris_update_ms_exp6 = [0 1 5 10 20 50 100];
        cfg.ITER_EXP7 = 20;
        cfg.exp7_random_trials = 50;
        % Fine block search is deliberately a high-compute reference:
        % it uses the same Stage-I seed and true objective, but a denser
        % derivative-free phase grid. This makes the performance/runtime
        % tradeoff interpretable rather than weakening the comparator.
        cfg.exp7_bcs_blocks = 8;
        cfg.exp7_bcs_phase_levels = 16;
        cfg.exp7_bcs_sweeps = 2;
        % Stage I uses Riemannian conjugate-gradient candidates and adaptive
        % proxy-screened exact evaluation. Every accepted point is checked by
        % the true objective.
        cfg.proposed_mechanism_fallback = true;
        cfg.proposed_stage1_geodesic_fractions = [1 0.5 0.25];
        cfg.proposed_stage1_exact_topk_initial = 2;
        cfg.proposed_stage1_exact_topk_max = 4;
        cfg.proposed_stage1_expand_on_stall = true;
        % Legacy polish settings retained for the fixed-RET/RIS sensing
        % baseline helper. They are not the final proposed V2 Stage-II.
        cfg.proposed_polish_blocks_schedule = [8 16];
        cfg.proposed_polish_phase_step_schedule = { ...
            [22.5 11.25],[11.25 5.625]};
        cfg.proposed_polish_active_blocks_schedule = [4 6];
        cfg.proposed_polish_exact_topk_initial_schedule = [2 2];
        cfg.proposed_polish_exact_topk_max_schedule = [4 4];
        cfg.proposed_polish_ret_half_window = 1;
        cfg.proposed_polish_ret_exact_topk_initial = 1;
        cfg.proposed_polish_ret_exact_topk_max = 2;
        cfg.proposed_polish_try_fine_on_stall = false;
        % Ignore numerically tiny accepted updates below 0.005 dB.
        cfg.proposed_polish_min_gain_dB = 5e-3;
        % MA-RCG-MF v2 Stage-II. Diagnosis of R3 showed that ~89% of the
        % MA-vs-Block gap came from restrictive candidate generation rather
        % than Top-K screening. V2 therefore visits every RIS block
        % sequentially, updates the incumbent immediately after each accepted
        % block move, and adapts a per-block trust radius. A full-circle
        % proxy-only coarse probe is used in the first sweep; only a small
        % adaptive Top-K subset receives the exact communication/QoS solve.
        cfg.proposed_stage2_variant = 'sequential_mf_atr';
        cfg.v2_blocks = 8;
        cfg.v2_max_sweeps = 3;
        cfg.v2_global_probe_levels = 16;
        cfg.v2_global_probe_sweeps = 1;
        cfg.v2_trust_init_deg = 45;
        cfg.v2_trust_min_deg = 2.8125;
        cfg.v2_trust_max_deg = 135;
        cfg.v2_trust_expand = 1.5;
        cfg.v2_trust_shrink = 0.5;
        cfg.v2_local_fractions = [1 0.5 0.25];
        cfg.v2_exact_topk_initial = 2;
        cfg.v2_exact_topk_max = 4;
        cfg.v2_expand_budget_on_stall = true;
        cfg.v2_proxy_tie_gap_dB = 0.05;
        cfg.v2_ret_half_window = 1;
        cfg.v2_ret_exact_topk_initial = 1;
        cfg.v2_ret_exact_topk_max = 2;
        % Accept physically meaningful improvements down to 0.002 dB.
        cfg.v2_min_gain_dB = 2e-3;
        cfg.v2_sweep_stop_gain_dB = 5e-3;
        cfg.v2_max_noaccept_sweeps = 2;
        % Engineering validation experiments.
        cfg.ITER_EXP8 = 30;
        cfg.ris_x_range_exp8 = [40 60 80 100 120];
        cfg.ITER_EXP9 = 100;
        cfg.exp9_bootstrap_reps = 2000;
        % Predeclared CSI stress grid. The -8/-5 dB points deliberately push
        % the fixed perfect-CSI cohort toward its outage boundary.
        cfg.csi_nmse_dB_range_exp9 = [-Inf -30 -20 -15 -10 -8 -5];
        cfg.csi_error_range_exp9 = 10.^(cfg.csi_nmse_dB_range_exp9/20);
        cfg.ITER_EXP10 = 50;
        cfg.ris_phase_bits_exp10 = [1 2 3 Inf];
        cfg.ITER_EXP11 = 30;
        cfg.ret_resolution_deg_exp11 = [5 2 1 0.5];
        cfg.ITER_EXP12 = 20;
    otherwise
        error('Unknown run mode: %s. Use ''fast'' or ''paper''.', run_mode);
end
cfg.num_workers = min(cfg.num_workers, cfg.ITER);

%% Common system size
cfg.M = 16;
cfg.N = 64;
cfg.K = 4;
cfg.T = 3;

%% Experiment axes
% Experiment 1 minimizes P_comm. This value is only a generous feasibility
% cap, not the quantity being optimized.
cfg.P_total_exp1 = 40; % W

cfg.P_total_dBm_range_exp2 = 24:2:32;
cfg.P_total_range_exp2 = 10.^((cfg.P_total_dBm_range_exp2 - 30)/10);
cfg.gamma_exp2_dB = 5;

cfg.P_total_exp3_dBm = 32;
cfg.P_total_exp3 = 10^((cfg.P_total_exp3_dBm - 30)/10);

cfg.radar_snr_target_dB_exp2 = [];
cfg.num_radar_snr_targets = 5;

%% Experiment 4: dual-domain ablation
% A moderate QoS point is selected so that all designs are normally feasible
% and the separate effects of power release and sensing-channel shaping can
% be observed without the high-QoS feasibility boundary dominating.
cfg.gamma_exp4_dB = 8;
cfg.P_total_exp4_dBm = 32;
cfg.P_total_exp4 = 10^((cfg.P_total_exp4_dBm-30)/10);

%% Experiment 5: scalability / complexity
cfg.gamma_exp5_dB = 8;
cfg.P_total_exp5_dBm = 32;
cfg.P_total_exp5 = 10^((cfg.P_total_exp5_dBm-30)/10);
% The performance curves are never post-processed to enforce monotonicity.
% Continuation and independent restarts only improve numerical robustness; the
% best candidate is selected using the exact MinPc or exact radar-SNR metric.
% If explicitly enabled for supplementary diagnostics, the timing benchmark
% is serial and uses the selected execution profile. CVX warm-up is discarded,
% and medians/IQRs are reported over repeated complete executions. The
% manuscript runtime comparison is taken from Exp.7/12, not this optional pass.

%% Experiment 6: multi-timescale mobility / stale-configuration loss
% RET is held over a slow traffic-control epoch. RIS is refreshed every
% tau_RIS milliseconds. Digital communication and sensing beamformers are
% re-optimized from current CSI at the fast timescale. The experiment
% isolates geometry-induced passive-configuration staleness.
cfg.gamma_exp6_dB = 8;
cfg.P_total_exp6_dBm = 32;
cfg.P_total_exp6 = 10^((cfg.P_total_exp6_dBm-30)/10);
cfg.exp6_user_motion_axis = 'x';
cfg.exp6_keep_RET_fixed = true;
cfg.exp6_fast_active_beam_update = true;
cfg.exp6_fresh_ris_trials = []; % optional cap on deterministic fresh-RIS candidates
cfg.exp6_zero_ms_calibration = true;
cfg.exp6_warm_start_from_stale = true;
% The stale phase is always retained as a feasible warm-start candidate.
% Additional current-snapshot candidates are deterministic target/communication
% focusing phases; unrelated random multi-starts are disabled.
cfg.exp6_warm_start_only = false; % optional diagnostic mode; false in all main experiments
cfg.exp6_clip_negative_loss = false;
cfg.exp6_numerical_loss_tol_dB = 1e-8;
cfg.exp6_bootstrap_seed = 260806;
cfg.exp6_report_displacement_curve = true;
% Mobility is applied to both communication vehicles and sensing targets.
% This is essential: moving only communication users leaves the target-side
% RIS phase almost unchanged and cannot reveal sensing-configuration aging.
cfg.exp6_move_sensing_targets = true;
cfg.exp6_target_speed_scale = 1.0;
cfg.exp6_target_motion_axis = 'x';
% Preserve the t=0 channel realization and apply the deterministic phase
% rotation caused by path-length changes. This yields a paired mobility test
% without inventing independent fading at every update period.
cfg.exp6_apply_relative_path_phase = true;
cfg.exp6_fresh_candidate_mode = 'warm_target_communication';
cfg.exp6_fresh_allow_random = false;
% Strict paired-mobility controls.  Every speed/update-period cell within an
% MC replicate is generated from one precomputed common-random-number state
% bank before any optimization is called.  No cell-specific random restart is
% permitted in the fresh-RIS benchmark.
cfg.exp6_pairing_version = 'strict-crn-state-bank-v1';
cfg.exp6_enforce_deterministic_fresh = true;
cfg.exp6_pairing_rebuild_check = true;
cfg.exp6_use_paired_bootstrap = true;
cfg.exp6_trend_audit = true;
cfg.exp6_resume_checkpoint = true;
cfg.exp6_checkpoint_batch = 10;
cfg.exp6_delete_checkpoint_on_success = true;
cfg.exp6_primary_metric = 'paired-endpoint-loss';
% Antithetic motion is available only as an optional sensitivity analysis; it
% is deliberately disabled for the manuscript so the MC unit remains one
% physical traffic realization rather than an averaged forward/reverse pair.
cfg.exp6_antithetic_motion = false;
% The main figure reports mean and tail radar-SNR aging loss. Communication-
% power change and a predeclared near-boundary outage stress point are retained
% only as supplementary diagnostics; their signs/values are not post-processed.
cfg.exp6_stress_gamma_dB = 12;
cfg.exp6_stress_P_total_dBm = 30;
cfg.exp6_stress_P_total = 10^((cfg.exp6_stress_P_total_dBm-30)/10);


%% Engineering-validation operating points
cfg.gamma_exp8_dB = 8;
cfg.P_total_exp8_dBm = 32;
cfg.P_total_exp8 = 10^((cfg.P_total_exp8_dBm-30)/10);
cfg.ris_y_exp8 = 8;

cfg.gamma_exp9_dB = 14;
cfg.P_total_exp9_dBm = 32;
cfg.P_total_exp9 = 10^((cfg.P_total_exp9_dBm-30)/10);
cfg.csi_error_model_exp9 = ['Hhat=sqrt(1-epsilon^2)H+epsilon E, with matched link-wise energy; ', ...
    'epsilon=10^(NMSE_target_dB/20), and Perfect CSI uses epsilon=0'];
cfg.csi_nmse_definition_exp9 = ['The plotted NMSE grid specifies the target error-power ratio epsilon^2. ', ...
    'Measured link-wise NMSE is also saved for verification.'];
cfg.exp9_revision = 'NMSE-grid-v3-fixed-cohort-predeclared-stress';
cfg.exp9_min_fixed_cohort_warning = 20;
cfg.csi_qos_tolerance_exp9 = 5e-4;

cfg.gamma_exp10_dB = 8;
cfg.P_total_exp10_dBm = 32;
cfg.P_total_exp10 = 10^((cfg.P_total_exp10_dBm-30)/10);

cfg.gamma_exp11_dB = 8;
cfg.P_total_exp11_dBm = 32;
cfg.P_total_exp11 = 10^((cfg.P_total_exp11_dBm-30)/10);
cfg.ret_local_search_span_deg_exp11 = 2;

cfg.gamma_exp12_dB = 8;
cfg.P_total_exp12_dBm = 32;
cfg.P_total_exp12 = 10^((cfg.P_total_exp12_dBm-30)/10);
cfg.significance_alpha_exp12 = 0.05;

%% CSI / target-information assumption used by the optimization
cfg.csi_assumption = 'perfect current active-link CSI and geometry-derived sensing CSI';
cfg.csi_user_links = {'BS-user','BS-RIS','RIS-user'};
cfg.target_information = 'target position/angle and large-scale reflection weights available from the preceding sensing epoch';
cfg.csi_future_work = 'robust design under channel-estimation, localization and target-state errors';

%% Carrier / wavelength (used by the mobility phase-aging test)
cfg.fc_Hz = 2.4e9;
cfg.c_light = 299792458;
cfg.lambda = cfg.c_light/cfg.fc_Hz;

%% Receiver-noise model
% The numerical noise level used in every experiment is derived rather than
% entered as an unexplained constant.  With -174 dBm/Hz thermal noise,
% 100 MHz effective bandwidth and 14 dB receiver noise figure, the resulting
% power is -80 dBm = 1e-11 W.  The sensing receiver uses the same effective
% bandwidth/noise figure for a transparent like-for-like SNR normalization.
cfg.noise_psd_dBm_Hz = -174;
cfg.comm_bandwidth_Hz = 100e6;
cfg.comm_noise_figure_dB = 14;
cfg.radar_bandwidth_Hz = 100e6;
cfg.radar_noise_figure_dB = 14;
cfg.sigma_c2_dBm = cfg.noise_psd_dBm_Hz + ...
    10*log10(cfg.comm_bandwidth_Hz) + cfg.comm_noise_figure_dB;
cfg.sigma_r2_dBm = cfg.noise_psd_dBm_Hz + ...
    10*log10(cfg.radar_bandwidth_Hz) + cfg.radar_noise_figure_dB;
cfg.sigma_c2 = 1e-3*10^(cfg.sigma_c2_dBm/10);
cfg.sigma_r2 = 1e-3*10^(cfg.sigma_r2_dBm/10);

%% Three-dimensional V2X geometry (metres)
cfg.bs_x = 0;  cfg.bs_y = 0;  cfg.h_BS = 25;
cfg.ris_x = 90; cfg.ris_y = 8; cfg.h_RIS = 5;
cfg.h_user = 1.5;
cfg.h_target = 1.5;

cfg.road_x_min = 60;
cfg.road_x_max = 150;
cfg.road_y_half = 4;

cfg.target_x = [70, 105, 140];
cfg.target_y = [-2, 1, 3];
% Equal effective target reflection/RCS is assumed in the main paper.
% target_rcs_linear may be replaced by measured or scenario-specific values;
% build_radar_matrix uses the normalized reflection weights below.
cfg.target_rcs_linear = ones(cfg.T,1);
cfg.target_weights = cfg.target_rcs_linear(:)/sum(cfg.target_rcs_linear);

%% Large-scale fading
cfg.C0_dB = -30;
cfg.C0 = 10^(cfg.C0_dB/10);
cfg.D0 = 1;
cfg.alpha_BR = 2.2;
cfg.alpha_Ru = 2.2;
cfg.alpha_Bu = 3.6;
cfg.alpha_Bt = 2.7;
cfg.alpha_Rt = 2.3;
cfg.beta_Ru_dB = 10;
cfg.beta_Ru = 10^(cfg.beta_Ru_dB/10);
cfg.direct_extra_loss_dB = 6;

%% RET electrical downtilt hardware assumption
% Sector-level common electronic downtilt: all active BS antenna ports share
% one common tilt variable. Per-user spatial separation is handled by the
% digital beamformers; no antenna element has an independent mechanical
% boresight variable as in per-antenna rotatable-antenna (RA) systems.
cfg.ret_hardware_type = 'sector-common electronic downtilt';
cfg.ret_control_granularity = 'one common angle for the full BS panel';
cfg.ret_angle_resolution_deg = 1;
cfg.ret_timescale = 'slow traffic/coverage control epoch';
cfg.theta_3dB = 65;
cfg.SLA_V = 30;
cfg.N_RET = 8;
cfg.d_RET_lambda = 0.5;

%% Fast exact communication-QoS inner solver
% 'hybrid' uses uplink-downlink duality/fixed-point first and automatically
% falls back to CVX/SOCP if the SINR/primal-dual validation fails. Set the
% RET_RIS_COMM_SOLVER environment variable to 'cvx' for legacy reproduction.
cfg.comm_solver_default = 'hybrid';
cfg.comm_duality_fixed_point_tol = 1e-7;
cfg.comm_duality_max_iter = 500;
cfg.comm_duality_primal_dual_gap_tol = 2e-4;
cfg.comm_solver_validation_power_gap_tol = 3e-3;
cfg.comm_solver_validation_sinr_tol = 5e-4;
cfg.comm_solver_runtime_accounting = 'complete executable path with one common backend policy';

%% AO / phase optimization
cfg.rho_init = 1;
cfg.c_rho = 0.8;
cfg.tol_obj = 1e-4;
cfg.tol_aux = 1e-4;
cfg.phase_max_iter = 100;
cfg.phase_grad_tol = 1e-6;

%% Fast-good proposed RET-RIS optimizer
% The exact structured sensing objective is J = J_comm + P_r*G_eff, with
% P_r=P_T-P_c^min and G_eff=lambda_max(Z'*C_r*Z). The implementation combines
% Riemannian conjugate-gradient (RCG) phase proposals, mechanism-aware proxy
% ranking, adaptive top-K exact evaluation, and V2 sequential all-block
% multi-fidelity refinement with per-block adaptive trust regions. The
% communication QoS subproblem is solved by a validated
% uplink-downlink-duality fixed-point method with automatic CVX fallback.
% Every passive-control update is accepted only after true-objective checking.
cfg.proposed_hybrid_comm_weight_min = 0.20;
cfg.proposed_hybrid_comm_weight_max = 0.80;
cfg.proposed_refine_accept_tol = 1e-10;
cfg.proposed_record_stage1 = true;
% Apply the same lightweight phase-only exact refinement to strong fixed-RET
% sensing baselines when they run in full optimization mode.  Mobility fresh-
% RIS updates remain restricted and therefore skip this baseline polish.
cfg.fixed_ris_enable_exact_polish = true;

%% Fixed power-split baselines for Experiment 2
cfg.alpha_fixed_list = [0.3, 0.5, 0.7];

%% Reproducibility
cfg.seed_base = 2025;
cfg.seed_stride = 1000;
% Common deterministic solver stream used by the shared Exp.2/Exp.3 sweep.
cfg.solver_point_seed_stride = 104729;
cfg.solver_seed_offset_sweep = 2500000;
% Exp.2/3 points are solved independently so an overlapping operating point
% is numerically identical in both experiments and is not path-dependent.
cfg.exp23_use_continuation = false;
cfg.ci_confidence = 0.95;

%% Derived horizontal geometry and azimuths
cfg.d_BR_horizontal = hypot(cfg.ris_x-cfg.bs_x, cfg.ris_y-cfg.bs_y);
cfg.bs2ris_az_deg = atan2d(cfg.ris_y-cfg.bs_y, cfg.ris_x-cfg.bs_x);
cfg.ris2bs_az_deg = atan2d(cfg.bs_y-cfg.ris_y, cfg.bs_x-cfg.ris_x);

cfg.d_Bt_horizontal = hypot(cfg.target_x-cfg.bs_x, cfg.target_y-cfg.bs_y);
cfg.d_Rt_horizontal = hypot(cfg.target_x-cfg.ris_x, cfg.target_y-cfg.ris_y);
cfg.target_az_deg = atan2d(cfg.target_y-cfg.bs_y, cfg.target_x-cfg.bs_x);
cfg.RIS2tar_deg = atan2d(cfg.target_y-cfg.ris_y, cfg.target_x-cfg.ris_x);

%% Three-dimensional propagation distances used in path loss
cfg.d_BR = sqrt(cfg.d_BR_horizontal.^2 + (cfg.h_BS-cfg.h_RIS).^2);
cfg.d_Bt_vec = sqrt(cfg.d_Bt_horizontal.^2 + (cfg.h_BS-cfg.h_target).^2);
cfg.d_Rt_vec = sqrt(cfg.d_Rt_horizontal.^2 + (cfg.h_RIS-cfg.h_target).^2);

%% Fair fixed-RET baseline
% The fixed baseline is not an arbitrary weak tilt. It is the sector-level
% geometry-centred downtilt pointing to the centre of the road service area,
% quantized to the same 1-deg hardware grid used by optimized RET.
cfg.fixed_ret_reference_x = 0.5*(cfg.road_x_min+cfg.road_x_max);
cfg.fixed_ret_reference_y = 0;
cfg.fixed_ret_raw_deg = atan2d(cfg.h_BS-cfg.h_user, ...
    hypot(cfg.fixed_ret_reference_x-cfg.bs_x, ...
          cfg.fixed_ret_reference_y-cfg.bs_y));
cfg.theta_ref = cfg.ret_angle_resolution_deg * ...
    round(cfg.fixed_ret_raw_deg/cfg.ret_angle_resolution_deg);
cfg.fixed_ret_selection = 'geometry-centred road-coverage angle, quantized to hardware grid';

%% Experiment 7: fair same-model algorithm comparison
cfg.gamma_exp7_dB = 8;
cfg.P_total_exp7_dBm = 32;
cfg.P_total_exp7 = 10^((cfg.P_total_exp7_dBm-30)/10);
cfg.exp7_compare_external_RA_code = false;
cfg.exp7_include_bcs_baseline = true;
cfg.exp7_bcs_name = 'Fine AO block-coordinate phase search';
cfg.exp7_external_RA_reason = ['Published RA-IRS work uses per-antenna 3-D ', ...
    'rotation, uplink sum-rate and a different channel/objective; no public ', ...
    'code was located. It is discussed in Related Work but not used as an ', ...
    'unfair numerical baseline.'];
end
