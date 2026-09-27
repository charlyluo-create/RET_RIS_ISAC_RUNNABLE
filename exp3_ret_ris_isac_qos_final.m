% =========================================================================
% Experiment 3 FINAL (Unified RET-RIS ISAC QoS Coupling Analysis)
%
% Purpose:
%   With a fixed total transmit-power budget, scan the communication SINR
%   requirement and quantify how communication QoS compresses the sensing
%   resources and radar performance.
%
% The script uses exactly the same M/N/K/T, road geometry, target positions,
% fading parameters, noise powers, RET model, RIS model and random seeds as
% the unified Experiments 1 and 2.
%
% Five schemes:
%   1) Fixed RET, No RIS
%   2) Fixed RET, RIS-MinPc
%   3) RET-MinPc, No RIS (MinPc-oriented)
%   4) RET-RIS-MinPc
%   5) RET-RIS-RadarSNR
%
% Optimization backbone:
%   - Wc: exact minimum-communication-power beamforming (duality FPI; CVX fallback)
%   - phi: penalty + Riemannian manifold optimization
%   - RET: complete grid for no-RIS initialization and local discrete search
%   - Wr: communication-null-space projection + principal eigenvector
% =========================================================================

clear; close all; clc;

%% ======================== Time And Progress Control =======================
% Independent timer handles prevent tic/toc calls inside CVX or Manopt from
% affecting the experiment-level timing results.
t_total = tic;
SHOW_PROGRESS = true;   % true: show progress, elapsed time and ETA

%% ======================== Shared Configuration ============================
env_mode=getenv('RET_RIS_RUN_MODE');
if isempty(env_mode); RUN_MODE='paper'; else; RUN_MODE=lower(env_mode); end
cfg = unified_ret_ris_isac_config(RUN_MODE);
FIG_POLICY=get_figure_policy(RUN_MODE);
SAVE_MAIN_FIGURE=FIG_POLICY.save_main;
SAVE_SUPP_FIGURES=FIG_POLICY.save_supplementary;
SAVE_FIGURES=FIG_POLICY.save_any;

ITER = cfg.ITER;
M = cfg.M; N = cfg.N; K = cfg.K; T = cfg.T;
gamma_dB_range = cfg.gamma_dB_range_exp3;
gamma_range = 10.^(gamma_dB_range/10);
numG = numel(gamma_dB_range);
P_total = cfg.P_total_exp3;
P_total_dBm = cfg.P_total_exp3_dBm;
num_workers = cfg.num_workers;

%% ======================== Path Setup =====================================
root = fileparts(mfilename('fullpath'));
if isempty(root); root = pwd; end
addpath(root);

% Create a new timestamped figure folder for every run so that previous
% figures are never overwritten.
run_timestamp = datestr(now,'yyyymmdd_HHMMSS');
figure_dir = fullfile(root, sprintf('Exp3_Figures_%s_%s', ...
    lower(RUN_MODE), run_timestamp));
if SAVE_FIGURES && ~exist(figure_dir,'dir')
    mkdir(figure_dir);
end

func_path = fullfile(root,'function');
manopt_path = fullfile(root,'manopt');
if exist(func_path,'dir'); addpath(func_path); end
if exist(manopt_path,'dir'); addpath(genpath(manopt_path)); end

required_functions = {'generate_channel','solve_comm_min_power', ...
    'opt_phi_min_pc','opt_phi_radar_gain','opt_a_min_pc', ...
    'build_radar_matrix','allocate_radar_power','build_effective_user_channel'};
for ii = 1:numel(required_functions)
    if exist(required_functions{ii},'file') ~= 2
        error('Required function not found: %s.m',required_functions{ii});
    end
end
if exist('cvx_begin','file') ~= 2
    error('CVX is required. Install CVX and run cvx_setup first.');
end

%% ======================== Parallel Configuration =========================
use_parallel = true;
if use_parallel
    try
        poolobj = gcp('nocreate');
        if ~isempty(poolobj) && poolobj.NumWorkers ~= num_workers
            delete(poolobj); poolobj = [];
        end
        if isempty(poolobj); poolobj = parpool(num_workers); end
        pctRunOnAll(['addpath(''' root ''');']);
        if exist(func_path,'dir'); pctRunOnAll(['addpath(''' func_path ''');']); end
        if exist(manopt_path,'dir')
            pctRunOnAll(['addpath(genpath(''' manopt_path '''));']);
        end
    catch ME
        warning('Parallel setup failed: %s. Falling back to serial.',ME.message);
        use_parallel = false;
    end
end

fprintf('\n=================================================================\n');
fprintf(' EXPERIMENT 3: UNIFIED RET-RIS ISAC QoS COUPLING\n');
fprintf(' Mode=%s, MC=%d, M=%d, N=%d, K=%d, T=%d\n', ...
    upper(RUN_MODE),ITER,M,N,K,T);
fprintf(' Fixed P_total=%.1f dBm (%.4f W), gamma=%g:%g:%g dB\n', ...
    P_total_dBm,P_total,gamma_dB_range(1), ...
    gamma_dB_range(2)-gamma_dB_range(1),gamma_dB_range(end));
fprintf(' Road targets: '); fprintf('(%.0f,%.0f) ',[cfg.target_x;cfg.target_y]); fprintf('m\n');
fprintf(' Same M/N/K/T, road, fading, RET, RIS and seeds as Experiments 1/2.\n');
if SAVE_FIGURES
    fprintf(' Figures will be saved to: %s\n',figure_dir);
end
fprintf('=================================================================\n\n');

%% ======================== Result Buffers =================================
S = 5;
Obj = nan(ITER,numG,S);          % weighted radar SNR, linear
Pcomm = nan(ITER,numG,S);        % minimum communication power, W
Pradar = nan(ITER,numG,S);       % residual radar power, W
Feasible = false(ITER,numG,S);
Theta = nan(ITER,numG,S);

%% ======================== Monte Carlo ====================================
fprintf('Starting Monte Carlo (%d iterations, %d QoS points)...\n',ITER,numG);

t_start = tic;

% Progress is counted at the granularity of one completed
% (Monte Carlo realization, QoS threshold) pair.
progressQueue = [];
if SHOW_PROGRESS
    try
        total_progress_steps = ITER * numG;
        progressQueue = parallel.pool.DataQueue;
        afterEach(progressQueue,@exp3_progress_bar);
        exp3_progress_bar(struct('action','reset','total',total_progress_steps));
    catch ME
        warning('Progress bar initialization failed: %s',ME.message);
        progressQueue = [];
    end
end

if use_parallel
    parfor iter = 1:ITER
        rng(cfg.seed_base + iter*cfg.seed_stride,'twister');
        seed=cfg.seed_base + iter*cfg.seed_stride;
        out = run_one_mc_qos(cfg,gamma_range,P_total,seed,func_path,progressQueue);
        Obj(iter,:,:) = reshape(out.Obj,[1,numG,S]);
        Pcomm(iter,:,:) = reshape(out.Pcomm,[1,numG,S]);
        Pradar(iter,:,:) = reshape(out.Pradar,[1,numG,S]);
        Feasible(iter,:,:) = reshape(out.Feasible,[1,numG,S]);
        Theta(iter,:,:) = reshape(out.Theta,[1,numG,S]);
    end
else
    for iter = 1:ITER
        rng(cfg.seed_base + iter*cfg.seed_stride,'twister');
        seed=cfg.seed_base + iter*cfg.seed_stride;
        out = run_one_mc_qos(cfg,gamma_range,P_total,seed,func_path,progressQueue);
        Obj(iter,:,:) = reshape(out.Obj,[1,numG,S]);
        Pcomm(iter,:,:) = reshape(out.Pcomm,[1,numG,S]);
        Pradar(iter,:,:) = reshape(out.Pradar,[1,numG,S]);
        Feasible(iter,:,:) = reshape(out.Feasible,[1,numG,S]);
        Theta(iter,:,:) = reshape(out.Theta,[1,numG,S]);
    end
end
runtime = toc(t_start);

%% ======================== Statistics =====================================
% Main performance curves are computed on realizations that are feasible for
% every scheme at the same QoS point. Scheme-specific feasibility remains a
% separate reliability metric.
CommonFeasible=all(Feasible,3);
CommonFeasibleRate=100*mean(CommonFeasible,1).';
Fea_rate=squeeze(mean(Feasible,1));

[Obj_dB,Obj_lo_dB,Obj_hi_dB,Obj_n,Obj_avg]= ...
    masked_mean_ci_to_db(Obj,CommonFeasible,1);
[Pcomm_dBm,Pcomm_lo_dBm,Pcomm_hi_dBm,Pcomm_n,Pcomm_avg]= ...
    masked_mean_ci_to_db(Pcomm,CommonFeasible,1e-3);
[Pradar_avg,Pradar_lo,Pradar_hi,Pradar_n]= ...
    masked_mean_ci_linear(Pradar,CommonFeasible);
[Theta_avg,Theta_lo,Theta_hi,Theta_n]= ...
    masked_mean_ci_linear(Theta,CommonFeasible);

Pcomm_ratio=Pcomm_avg/P_total;
Pradar_ratio=Pradar_avg/P_total;
RadarGain_dB=nan(numG,S); RadarGain_lo=nan(numG,S); RadarGain_hi=nan(numG,S);
RadarGain_n=zeros(numG,S); RadarGain_dB(:,1)=0;
for ss=2:S
    [RadarGain_dB(:,ss),RadarGain_lo(:,ss),RadarGain_hi(:,ss),RadarGain_n(:,ss)]= ...
        paired_db_statistics(Obj(:,:,ss),Obj(:,:,1), ...
        Feasible(:,:,ss)&Feasible(:,:,1));
end
labels = {'Fixed RET, No RIS', ...
          'Fixed RET, RIS-MinPc', ...
          'RET-MinPc, No RIS', ...
          'RET-RIS-MinPc', ...
          'RET-RIS-RadarSNR'};
markers = {'-o','-^','-d','-s','-p'};

if SAVE_MAIN_FIGURE
%% ======================== Main Paper Composite ===========================
% The manuscript uses one 2x2 QoS-coupling figure. Residual-power and RET-
% angle panels remain supplementary diagnostics. No curve smoothing or
% monotonic envelope is applied.
fig0 = figure('Position',[30 30 1450 900],'Color','w');
tl0 = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

nexttile;
for s=1:S
    errorbar(gamma_dB_range,Pcomm_dBm(:,s),Pcomm_dBm(:,s)-Pcomm_lo_dBm(:,s), ...
        Pcomm_hi_dBm(:,s)-Pcomm_dBm(:,s),markers{s},'LineWidth',1.5,'MarkerSize',7); hold on;
end
yline(P_total_dBm,'--k','P_T','HandleVisibility','off');
grid on; xlabel('Communication QoS threshold \gamma (dB)');
ylabel('Communication power P_c^{min} (dBm)');
title(sprintf('(a) Common-feasible communication pressure (min n=%d)',min(Pcomm_n(:)))); legend(labels,'Location','best');

nexttile;
for s=1:S
    errorbar(gamma_dB_range,Obj_dB(:,s),Obj_dB(:,s)-Obj_lo_dB(:,s), ...
        Obj_hi_dB(:,s)-Obj_dB(:,s),markers{s},'LineWidth',1.5,'MarkerSize',7); hold on;
end
grid on; xlabel('Communication QoS threshold \gamma (dB)');
ylabel('Weighted radar SNR (dB)');
title('(b) Common-feasible sensing degradation under QoS coupling');

nexttile;
for s=1:S
    plot(gamma_dB_range,100*Fea_rate(:,s),markers{s},'LineWidth',1.5,'MarkerSize',7); hold on;
end
plot(gamma_dB_range,CommonFeasibleRate,'--k','LineWidth',1.5);
grid on; ylim([0 105]); xlabel('Communication QoS threshold \gamma (dB)');
ylabel('Feasibility rate (%)'); title('(c) Scheme-specific and common feasibility');
legend([labels,{'Common feasible samples'}],'Location','best');

nexttile;
for s=1:S
    good=isfinite(Pcomm_ratio(:,s)) & isfinite(Obj_dB(:,s));
    plot(100*Pcomm_ratio(good,s),Obj_dB(good,s),markers{s}, ...
        'LineWidth',1.5,'MarkerSize',7); hold on;
end
grid on; xlabel('Communication-power share 100P_c/P_T (%)');
ylabel('Weighted radar SNR (dB)'); title('(d) Communication-sensing resource tradeoff');

sgtitle(tl0,sprintf('QoS-Coupled RET-RIS ISAC (P_T=%.1f dBm)',P_total_dBm));
if SAVE_MAIN_FIGURE
    save_exp3_figure(fig0,figure_dir,'00_Main_QoS_Coupling_Composite');
end

end

if SAVE_SUPP_FIGURES
%% ======================== Figure 1: Communication Pressure ================
fig1 = figure('Position',[60 60 920 610],'Color','w');
for s=1:S
    plot(gamma_dB_range,Pcomm_dBm(:,s),markers{s}, ...
        'LineWidth',1.6,'MarkerSize',8); hold on;
end
yline(P_total_dBm,'--k','P_T','LineWidth',1.2);
xlabel('Communication QoS threshold \gamma (dB)');
ylabel('Minimum communication power P_c^{min} (dBm)');
legend([labels,{'Total-power budget'}],'Location','northwest');
grid on; title('Communication Power Pressure Under Increasing QoS');
if SAVE_SUPP_FIGURES
    save_exp3_figure(fig1,figure_dir,'S1_Communication_Power_Pressure');
end

%% ======================== Figure 2: Radar Resource Share ==================
fig2 = figure('Position',[80 80 920 610],'Color','w');
for s=1:S
    plot(gamma_dB_range,100*Pradar_ratio(:,s),markers{s}, ...
        'LineWidth',1.6,'MarkerSize',8); hold on;
end
xlabel('Communication QoS threshold \gamma (dB)');
ylabel('Residual radar-power share 100P_r/P_T (%)');
legend(labels,'Location','southwest'); grid on;
title('Communication QoS Compresses the Radar Power Budget');
if SAVE_SUPP_FIGURES
    save_exp3_figure(fig2,figure_dir,'S2_Residual_Radar_Power_Share');
end

%% ======================== Figure 3: Radar SNR =============================
fig3 = figure('Position',[100 100 920 610],'Color','w');
for s=1:S
    plot(gamma_dB_range,Obj_dB(:,s),markers{s}, ...
        'LineWidth',1.6,'MarkerSize',8); hold on;
end
xlabel('Communication QoS threshold \gamma (dB)');
ylabel('Weighted radar SNR (dB)');
legend(labels,'Location','southwest'); grid on;
title('Radar Performance Under Communication QoS Coupling');
if SAVE_SUPP_FIGURES
    save_exp3_figure(fig3,figure_dir,'S3_Weighted_Radar_SNR');
end

%% ======================== Figure 4: Feasibility Rate ======================
fig4 = figure('Position',[120 120 920 610],'Color','w');
for s=1:S
    plot(gamma_dB_range,100*Fea_rate(:,s),markers{s}, ...
        'LineWidth',1.6,'MarkerSize',8); hold on;
end
xlabel('Communication QoS threshold \gamma (dB)');
ylabel('Feasibility rate (%)'); ylim([0 105]);
legend(labels,'Location','southwest'); grid on;
title(sprintf('System Feasibility at P_T = %.1f dBm',P_total_dBm));
if SAVE_SUPP_FIGURES
    save_exp3_figure(fig4,figure_dir,'S4_System_Feasibility_Rate');
end

%% ======================== Figure 5: Communication-Sensing Tradeoff ========
fig5 = figure('Position',[140 140 920 610],'Color','w');
for s=1:S
    good = isfinite(Pcomm_ratio(:,s)) & isfinite(Obj_dB(:,s));
    plot(100*Pcomm_ratio(good,s),Obj_dB(good,s),markers{s}, ...
        'LineWidth',1.6,'MarkerSize',8); hold on;
end
xlabel('Communication-power share 100P_c/P_T (%)');
ylabel('Weighted radar SNR (dB)');
legend(labels,'Location','southwest'); grid on;
title('Communication-Sensing Resource Tradeoff');
if SAVE_SUPP_FIGURES
    save_exp3_figure(fig5,figure_dir,'S5_Communication_Sensing_Tradeoff');
end

%% ======================== Figure 6: Optimal RET Tilt ======================
fig6 = figure('Position',[160 160 1000 450],'Color','w');
all_theta_plot=Theta_avg(:,3:5); all_theta_plot=all_theta_plot(isfinite(all_theta_plot));
if isempty(all_theta_plot)
    theta_ylim=[0 35];
else
    theta_ylim=[floor(min(all_theta_plot)-0.5),ceil(max(all_theta_plot)+0.5)];
    if diff(theta_ylim)<3; theta_ylim=mean(theta_ylim)+[-1.5 1.5]; end
end
subplot(1,3,1);
plot(gamma_dB_range,Theta_avg(:,3),'-d','LineWidth',1.5,'MarkerSize',7); grid on;
xlabel('\gamma (dB)'); ylabel('Optimal RET tilt (deg)'); title('RET-NoRIS-MinPc'); ylim(theta_ylim);
subplot(1,3,2);
plot(gamma_dB_range,Theta_avg(:,4),'-s','LineWidth',1.5,'MarkerSize',7); grid on;
xlabel('\gamma (dB)'); ylabel('Optimal RET tilt (deg)'); title('RET-RIS-MinPc'); ylim(theta_ylim);
subplot(1,3,3);
plot(gamma_dB_range,Theta_avg(:,5),'-p','LineWidth',1.5,'MarkerSize',7); grid on;
xlabel('\gamma (dB)'); ylabel('Optimal RET tilt (deg)'); title('RET-RIS-RadarSNR'); ylim(theta_ylim);
if SAVE_SUPP_FIGURES
    save_exp3_figure(fig6,figure_dir,'S6_Optimal_RET_Tilt_Angles');
end


%% ======================== Figure 7: Paired SNR Gains =====================
fig7=figure('Position',[180 180 950 620],'Color','w'); hold on;
for ss=2:S
    errorbar(gamma_dB_range,RadarGain_dB(:,ss), ...
        RadarGain_dB(:,ss)-RadarGain_lo(:,ss), ...
        RadarGain_hi(:,ss)-RadarGain_dB(:,ss),markers{ss}, ...
        'LineWidth',1.5,'MarkerSize',7);
end
yline(0,'--k','HandleVisibility','off'); grid on;
xlabel('Communication QoS threshold \gamma (dB)');
ylabel('Paired radar-SNR gain over fixed RET/no RIS (dB)');
legend(labels(2:S),'Location','best');
title('Realization-Paired Sensing Gain under Increasing QoS');
if SAVE_SUPP_FIGURES
    save_exp3_figure(fig7,figure_dir,'S7_Paired_Radar_SNR_Gain');
end

end

%% ======================== Tables / Save ==================================
fprintf('\nAverage minimum communication power (dBm):\n');
disp(array2table([gamma_dB_range(:),Pcomm_dBm], ...
    'VariableNames',{'gamma_dB','Fixed_NoRIS','Fixed_RIS_MinPc', ...
    'RET_NoRIS','RET_RIS_MinPc','RET_RIS_RadarSNR'}));

fprintf('\nAverage residual radar-power share:\n');
disp(array2table([gamma_dB_range(:),Pradar_ratio], ...
    'VariableNames',{'gamma_dB','Fixed_NoRIS','Fixed_RIS_MinPc', ...
    'RET_NoRIS','RET_RIS_MinPc','RET_RIS_RadarSNR'}));

fprintf('\nAverage weighted radar SNR (dB):\n');
disp(array2table([gamma_dB_range(:),Obj_dB], ...
    'VariableNames',{'gamma_dB','Fixed_NoRIS','Fixed_RIS_MinPc', ...
    'RET_NoRIS','RET_RIS_MinPc','RET_RIS_RadarSNR'}));

fprintf('\nRadar-SNR gain relative to fixed/no-RIS (dB):\n');
disp(array2table([gamma_dB_range(:),RadarGain_dB(:,2:5)], ...
    'VariableNames',{'gamma_dB','Fixed_RIS_MinPc','RET_NoRIS', ...
    'RET_RIS_MinPc','RET_RIS_RadarSNR'}));

fprintf('\nCommon-feasible rate across all schemes (%%):\n');
disp(array2table([gamma_dB_range(:),CommonFeasibleRate], ...
    'VariableNames',{'gamma_dB','CommonFeasible_percent'}));

fprintf('\nFeasibility rate:\n');
disp(array2table([gamma_dB_range(:),100*Fea_rate], ...
    'VariableNames',{'gamma_dB','Fixed_NoRIS','Fixed_RIS_MinPc', ...
    'RET_NoRIS','RET_RIS_MinPc','RET_RIS_RadarSNR'}));

results_dir = fullfile(root,'results');
if ~exist(results_dir,'dir'); mkdir(results_dir); end
save(fullfile(results_dir,sprintf('exp3_qos_revised_%s.mat',lower(RUN_MODE))), ...
    'cfg','gamma_dB_range','gamma_range','P_total','P_total_dBm', ...
    'Obj','Pcomm','Pradar','Feasible','Theta','CommonFeasible', ...
    'CommonFeasibleRate','Obj_avg','Obj_dB','Obj_lo_dB','Obj_hi_dB','Obj_n', ...
    'Pcomm_avg','Pcomm_dBm','Pcomm_lo_dBm','Pcomm_hi_dBm','Pcomm_n', ...
    'Pradar_avg','Pradar_ratio','Fea_rate','Theta_avg','Theta_lo','Theta_hi', ...
    'RadarGain_dB','RadarGain_lo','RadarGain_hi','RadarGain_n');
Tmain=array2table([gamma_dB_range(:),CommonFeasibleRate,Obj_dB,Pcomm_dBm,100*Fea_rate], ...
    'VariableNames',{'gamma_dB','CommonFeasible_percent', ...
    'Obj1_dB','Obj2_dB','Obj3_dB','Obj4_dB','Obj5_dB', ...
    'Pc1_dBm','Pc2_dBm','Pc3_dBm','Pc4_dBm','Pc5_dBm', ...
    'Fea1_percent','Fea2_percent','Fea3_percent','Fea4_percent','Fea5_percent'});
Tgain=array2table([gamma_dB_range(:),RadarGain_dB(:,2:S), ...
    RadarGain_lo(:,2:S),RadarGain_hi(:,2:S),RadarGain_n(:,2:S)], ...
    'VariableNames',{'gamma_dB', ...
    'Gain2_dB','Gain3_dB','Gain4_dB','Gain5_dB', ...
    'Gain2_Low','Gain3_Low','Gain4_Low','Gain5_Low', ...
    'Gain2_High','Gain3_High','Gain4_High','Gain5_High', ...
    'Gain2_N','Gain3_N','Gain4_N','Gain5_N'});
writetable(Tmain,fullfile(results_dir,sprintf('exp3_main_statistics_%s.csv',lower(RUN_MODE))));
writetable(Tgain,fullfile(results_dir,sprintf('exp3_paired_gains_%s.csv',lower(RUN_MODE))));
total_time = toc(t_total);
fprintf('\n');
fprintf('==================================================\n');
fprintf('  SIMULATION COMPLETED (Exp3 - RET-RIS QoS)\n');
fprintf('  ------------------------------------------------\n');
fprintf('  Monte Carlo time:  %8.1f s  (%6.2f min / %5.2f h)\n', ...
    runtime,runtime/60,runtime/3600);
fprintf('  Total runtime:     %8.1f s  (%6.2f min / %5.2f h)\n', ...
    total_time,total_time/60,total_time/3600);
fprintf('  Avg per MC run:    %8.1f s  (%6.2f min)\n', ...
    runtime/ITER,runtime/ITER/60);
fprintf('  Avg per QoS task:  %8.1f s\n',runtime/(ITER*numG));
fprintf('  Mode:              %s\n',upper(RUN_MODE));
if use_parallel
    fprintf('  Workers:           %d\n',num_workers);
end
if SAVE_FIGURES
    fprintf('  Figure folder:     %s\n',figure_dir);
end
fprintf('==================================================\n');

%% ========================================================================
%% Local Functions
%% ========================================================================

function out = run_one_mc_qos(cfg,gamma_range,P_total,seed,func_path,progressQueue)
% One paired roadside realization.  Exp. 3 shares the exact same five-scheme
% implementation and warm-start policy as Exp. 2.
if ~isempty(func_path) && exist(func_path,'dir'); addpath(func_path); end
scen=generate_isac_scenario(cfg,cfg.N,cfg.K,seed,[],[],struct());
solver_seed=seed+cfg.solver_seed_offset_sweep;
sw=evaluate_ret_ris_sweep(cfg,scen,gamma_range,P_total,solver_seed, ...
    struct('use_continuation',cfg.exp23_use_continuation, ...
    'seed_policy','operating_point'));
out.Obj=sw.Obj; out.Pcomm=sw.Pcomm; out.Pradar=sw.Pradar;
out.Feasible=sw.Feasible; out.Theta=sw.Theta;
for gIdx=1:numel(gamma_range)
    if ~isempty(progressQueue); send(progressQueue,1); end
end
end

function y=lin2db_safe(x)
y=nan(size(x)); valid=isfinite(x)&x>0; y(valid)=10*log10(x(valid));
end

function y=watt2dbm(x)
y=nan(size(x)); valid=isfinite(x)&x>0; y(valid)=10*log10(x(valid)/1e-3);
end

function m=nanmean_cols(X)
m=nan(1,size(X,2));
for ii=1:size(X,2)
    v=isfinite(X(:,ii)); if any(v); m(ii)=mean(X(v,ii)); end
end
end

%% ======================== Figure Save Helper ==============================
function save_exp3_figure(fig_handle,figure_dir,base_name)
% SAVE_EXP3_FIGURE Save one figure as a 300-dpi PNG and editable MATLAB FIG.

if ~exist(figure_dir,'dir')
    mkdir(figure_dir);
end

set(fig_handle,'Color','w');
drawnow;

png_file = fullfile(figure_dir,[base_name '.png']);
fig_file = fullfile(figure_dir,[base_name '.fig']);

try
    exportgraphics(fig_handle,png_file,'Resolution',300);
catch
    print(fig_handle,png_file,'-dpng','-r300');
end

try
    savefig(fig_handle,fig_file);
catch ME
    warning('Could not save editable FIG file %s: %s',fig_file,ME.message);
end

fprintf('Saved figure: %s\n',png_file);
end

%% ======================== Progress And Time Helpers =======================
function exp3_progress_bar(message)
% EXP3_PROGRESS_BAR Command-window progress bar with elapsed time and ETA.
% The callback executes on the MATLAB client. Workers only send scalar 1.

persistent total_steps completed_steps start_timer last_line_length;

if isstruct(message)
    if isfield(message,'action') && strcmp(message.action,'reset')
        total_steps = max(1,double(message.total));
        completed_steps = 0;
        start_timer = tic;
        last_line_length = 0;
        print_progress_line(false);
    end
    return;
end

if isempty(total_steps) || total_steps <= 0
    return;
end

completed_steps = min(total_steps, ...
    completed_steps + max(0,double(message)));
print_progress_line(completed_steps >= total_steps);

if completed_steps >= total_steps
    total_steps = [];
    completed_steps = [];
    start_timer = [];
    last_line_length = [];
end

    function print_progress_line(is_finished)
        elapsed_sec = toc(start_timer);

        if completed_steps > 0
            eta_sec = elapsed_sec / completed_steps * ...
                (total_steps - completed_steps);
        else
            eta_sec = NaN;
        end

        percentage = 100 * completed_steps / total_steps;
        bar_length = 30;
        filled_length = round(bar_length * completed_steps / total_steps);
        filled_length = min(max(filled_length,0),bar_length);
        bar_text = [repmat('=',1,filled_length), ...
                    repmat(' ',1,bar_length-filled_length)];

        elapsed_text = format_exp3_time(elapsed_sec);
        eta_text = format_exp3_time(eta_sec);

        line_text = sprintf(['Progress: [%s] %d/%d (%5.1f%%)  ', ...
            'Elapsed %s  ETA %s'], ...
            bar_text,completed_steps,total_steps,percentage, ...
            elapsed_text,eta_text);

        if isempty(last_line_length)
            last_line_length = 0;
        end
        if length(line_text) < last_line_length
            line_text = [line_text, ...
                repmat(' ',1,last_line_length-length(line_text))];
        end
        last_line_length = length(line_text);

        fprintf('\r%s',line_text);
        drawnow limitrate;

        if is_finished
            fprintf('\n');
        end
    end
end

function time_text = format_exp3_time(time_seconds)
% FORMAT_EXP3_TIME Convert seconds to HH:MM:SS for progress display.

if ~isfinite(time_seconds) || time_seconds < 0
    time_text = '--:--:--';
    return;
end

hours_value = floor(time_seconds/3600);
minutes_value = floor(mod(time_seconds,3600)/60);
seconds_value = floor(mod(time_seconds,60));

time_text = sprintf('%02d:%02d:%02d', ...
    hours_value,minutes_value,seconds_value);
end

