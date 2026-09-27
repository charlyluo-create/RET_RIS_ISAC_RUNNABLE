function run_all_paper_experiments(start_exp)
%RUN_ALL_PAPER_EXPERIMENTS Minimal clean entry point for the PAPER run.
%
% Normal use:
%   run_all_paper_experiments
%
% Resume after an interrupted run (earlier result MAT files must still exist):
%   run_all_paper_experiments(6)
%
% This clean package deliberately contains only the experiment pipeline,
% figure/table exporters, shared configuration, and the transitive helper
% functions required by those files. Legacy validators, patch scripts,
% one-off runners, manifests, and duplicate implementations are omitted.

if nargin < 1 || isempty(start_exp)
    start_exp = 1;
end
validateattributes(start_exp, {'numeric'}, ...
    {'scalar','integer','>=',1,'<=',12}, mfilename, 'start_exp');

root = fileparts(mfilename('fullpath'));
if isempty(root); root = pwd; end
func_dir = fullfile(root,'function');

% Put this release first on the MATLAB path so older copies elsewhere do not
% shadow the files in this folder.
addpath(root,'-begin');
addpath(func_dir,'-begin');
rehash path;

minimal_preflight(root, func_dir);

old_mode    = getenv('RET_RIS_RUN_MODE');
old_profile = getenv('RET_RIS_FIGURE_PROFILE');
old_exp5    = getenv('RET_RIS_EXP5_PROFILE');
old_skip    = getenv('RET_RIS_EXP5_SKIP_TIMING');
old_force   = getenv('RET_RIS_EXP5_FORCE_TIMING');
old_solver  = getenv('RET_RIS_COMM_SOLVER');
old_visible = get(groot,'defaultFigureVisible');
cleanupObj = onCleanup(@()restore_state(old_mode,old_profile,old_exp5, ...
    old_skip,old_force,old_solver,old_visible)); %#ok<NASGU>

setenv('RET_RIS_RUN_MODE','paper');
setenv('RET_RIS_FIGURE_PROFILE','paper_master');
setenv('RET_RIS_EXP5_PROFILE','balanced');
setenv('RET_RIS_EXP5_SKIP_TIMING','1');
setenv('RET_RIS_EXP5_FORCE_TIMING','0');
setenv('RET_RIS_COMM_SOLVER','hybrid');
set(groot,'defaultFigureVisible','off');

cfg = unified_ret_ris_isac_config('paper');
fprintf('\n============================================================\n');
fprintf('RET-RIS ISAC CLEAN PAPER PIPELINE\n');
fprintf('Release : %s\n', cfg.paper_release_version);
fprintf('Workers : %d (edit cfg.num_workers only if your machine supports more)\n', cfg.num_workers);
fprintf('Noise   : %.1f dBm communication, %.1f dBm sensing\n', ...
    10*log10(cfg.sigma_c2/1e-3), 10*log10(cfg.sigma_r2/1e-3));
fprintf('Exp.6   : %d MC, %d paired-bootstrap repetitions\n', ...
    cfg.ITER_EXP6, cfg.exp6_bootstrap_reps);
fprintf('Start   : experiment %d of 12\n', start_exp);
fprintf('============================================================\n');

fprintf('\n===== COMMUNICATION SOLVER DUALITY-vs-CVX VALIDATION =====\n');
validate_comm_solver_duality();

files = { ...
    'exp1_ret_ris_communication_final.m', ...
    'exp2_ret_ris_isac_power_final.m', ...
    'exp3_ret_ris_isac_qos_final.m', ...
    'exp4_ret_ris_ablation_final.m', ...
    'exp5_ret_ris_scalability_final.m', ...
    'exp6_ret_ris_mobility_multitimescale_final.m', ...
    'exp7_algorithm_comparison_convergence_final.m', ...
    'exp8_ris_deployment_final.m', ...
    'exp9_csi_uncertainty_final.m', ...
    'exp10_ris_phase_quantization_final.m', ...
    'exp11_ret_resolution_final.m', ...
    'exp12_complexity_significance_final.m'};

for ii = start_exp:numel(files)
    fprintf('\n===== PAPER %d/%d: %s =====\n', ii, numel(files), files{ii});
    run_isolated(fullfile(root,files{ii}));
    close all force;
end

fprintf('\n===== MA-RCG-MF V2 POST-RUN AUDIT =====\n');
run_isolated(fullfile(root,'exp7_v2_postrun_audit.m'));

fprintf('\n===== EXPORTING FINAL TABLES =====\n');
run_isolated(fullfile(root,'export_all_paper_statistics.m'));

fprintf('\n===== BUILDING SIX MAIN FIGURES =====\n');
figure_dir = generate_paper_figures('paper');
close all force;

fprintf('\n============================================================\n');
fprintf('PAPER pipeline completed successfully.\n');
fprintf('Results : %s\n', fullfile(root,'results'));
fprintf('Figures : %s\n', figure_dir);
fprintf('Main    : 6 composite figures / 22 panels\n');
fprintf('============================================================\n');
end

function minimal_preflight(root, func_dir)
% Embedded preflight: avoids a separate family of validator/patch files.
if ~exist(func_dir,'dir')
    error('Missing helper folder: %s',func_dir);
end
if exist('cvx_begin','file') ~= 2
    error(['CVX is not on the MATLAB path. Run cvx_setup first, then rerun ' ...
        'run_all_paper_experiments.']);
end

required_top = { ...
    'unified_ret_ris_isac_config.m','get_figure_policy.m', ...
    'generate_paper_figures.m','export_all_paper_statistics.m', ...
    'validate_comm_solver_duality.m', ...
    'exp7_v2_postrun_audit.m', ...
    'exp1_ret_ris_communication_final.m','exp2_ret_ris_isac_power_final.m', ...
    'exp3_ret_ris_isac_qos_final.m','exp4_ret_ris_ablation_final.m', ...
    'exp5_ret_ris_scalability_final.m','exp6_ret_ris_mobility_multitimescale_final.m', ...
    'exp7_algorithm_comparison_convergence_final.m','exp8_ris_deployment_final.m', ...
    'exp9_csi_uncertainty_final.m','exp10_ris_phase_quantization_final.m', ...
    'exp11_ret_resolution_final.m','exp12_complexity_significance_final.m'};
for ii=1:numel(required_top)
    p=fullfile(root,required_top{ii});
    if exist(p,'file')~=2
        error('Clean release is incomplete. Missing: %s',p);
    end
end

required_fun = { ...
    'allocate_radar_power','apply_relative_motion_phase','bootstrap_mean_ci_to_db', ...
    'bootstrap_quantile_grid','bootstrap_quantile_grid_paired','build_effective_user_channel', ...
    'build_radar_matrix','comm_solver_stats_add','comm_solver_stats_delta', ...
    'comm_solver_stats_summary','comm_solver_telemetry','empirical_quantile', ...
    'evaluate_external_beamformers','evaluate_fixed_ret_ris_config', ...
    'evaluate_fixed_split_configuration','evaluate_radar_only_reference', ...
    'evaluate_ret_ris_schemes','evaluate_ret_ris_sweep','generate_channel', ...
    'generate_isac_scenario','masked_mean_ci_linear','masked_mean_ci_to_db', ...
    'opt_a_min_pc','opt_phi_min_pc','opt_phi_radar_gain', ...
    'optimize_fixed_ret_ris_radar','optimize_ret_only_radar', ...
    'optimize_ret_ris_block_coordinate','paired_adjacent_trend_stats', ...
    'paired_db_statistics','paired_significance_test','perturb_scenario_csi', ...
    'quantize_ris_phase','refine_fixed_ret_ris_exact','refresh_ret_ris_geometry', ...
    'ret_ris_operating_point_seed','ret_vertical_gain','solve_comm_min_power', ...
    'solve_comm_min_power_cvx_backend','solve_comm_min_power_duality', ...
    'student_t_critical','unit_modulus_rcg'};
for ii=1:numel(required_fun)
    if exist(required_fun{ii},'file')~=2
        error('Missing required helper function: %s.m',required_fun{ii});
    end
end

cfg = unified_ret_ris_isac_config('paper');
assert(strcmp(cfg.paper_release_version,'2026-08-08-ma-rcg-mf-v2-sequential-mf-atr'), ...
    'Unexpected configuration release identifier.');
assert(isfield(cfg,'exp6_pairing_revision') && ...
    strcmp(cfg.exp6_pairing_revision,'2026-08-08-exp6-strict-paired-crn-v1'), ...
    'Strict-paired Exp.6 configuration is missing.');
assert(strcmpi(cfg.proposed_stage2_variant,'sequential_mf_atr'), ...
    'MA-RCG-MF V2 Stage-II configuration is missing.');
assert(any(cfg.K_range_exp5==10) && isequal(cfg.ris_phase_bits_exp10,[1 2 3 Inf]), ...
    'Experiment grids are inconsistent.');

% Detect, but do not fail on, duplicate copies elsewhere. This folder was
% already put first on the path, so the clean copy wins deterministically.
for name = {'unified_ret_ris_isac_config','generate_paper_figures'}
    allhits = which(name{1},'-all');
    if ischar(allhits); allhits={allhits}; end
    if numel(allhits)>1
        warning('Multiple copies of %s are on the MATLAB path. The clean release copy is used first.',name{1});
    end
end
fprintf('Minimal clean preflight passed.\n');
end

function run_isolated(file)
% Run scripts in the base workspace so their internal CLEAR statements do not
% destroy variables belonging to this driver function.
cmd = sprintf('run(''%s'');',strrep(file,'''',''''''));
evalin('base',cmd);
end

function restore_state(old_mode,old_profile,old_exp5,old_skip,old_force,old_solver,old_visible)
setenv('RET_RIS_RUN_MODE',old_mode);
setenv('RET_RIS_FIGURE_PROFILE',old_profile);
setenv('RET_RIS_EXP5_PROFILE',old_exp5);
setenv('RET_RIS_EXP5_SKIP_TIMING',old_skip);
setenv('RET_RIS_EXP5_FORCE_TIMING',old_force);
setenv('RET_RIS_COMM_SOLVER',old_solver);
set(groot,'defaultFigureVisible',old_visible);
end
