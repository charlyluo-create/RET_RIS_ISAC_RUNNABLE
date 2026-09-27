% =========================================================================
% Experiment 2 FINAL (Unified MinPc and Radar-SNR Comparison)
% RET-RIS-aided ISAC radar power reallocation and equal-sensing comparison
%
% Improvements:
%   1) total-power range is comparable with P_comm (24--32 dBm);
%   2) sensing targets lie on the same V2X road as communication users;
%   3) RIS multi-start includes target- and communication-focused phases;
%   4) RIS/RET updates are accepted by exact four-path radar SNR;
%   5) additional radar power is plotted in mW;
%   6) minimum total power at equal radar SNR is reported.
% =========================================================================

clear; close all; clc;

% ======================== Time And Progress Control ========================
t_total = tic;          % Whole-program timer (independent of Manopt/CVX tic/toc)
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
range_P = cfg.P_total_range_exp2;
range_P_dBm = cfg.P_total_dBm_range_exp2;
numP = numel(range_P);
gamma_dB = cfg.gamma_exp2_dB;
gamma = 10^(gamma_dB/10);
numAlpha = numel(cfg.alpha_fixed_list);
num_workers = cfg.num_workers;

%% ======================== Path Setup =====================================
root = fileparts(mfilename('fullpath'));
if isempty(root); root = pwd; end

%% ======================== Automatic Figure Folder =======================
% A new folder is created beside this script for every run.  The timestamp
% prevents FAST/PAPER runs from overwriting figures generated previously.
run_stamp = datestr(now, 'yyyymmdd_HHMMSS');
figure_dir = fullfile(root, sprintf('Exp2_Figures_%s_%s', ...
    lower(RUN_MODE), run_stamp));
if SAVE_FIGURES && ~exist(figure_dir, 'dir')
    mkdir(figure_dir);
end
if SAVE_FIGURES
    fprintf('Experiment 2 figures will be saved to:\n%s\n', figure_dir);
end

addpath(root);
func_path = fullfile(root,'function');
manopt_path = fullfile(root,'manopt');
if exist(func_path,'dir'); addpath(func_path); end
if exist(manopt_path,'dir'); addpath(genpath(manopt_path)); end
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
        if exist(manopt_path,'dir'); pctRunOnAll(['addpath(genpath(''' manopt_path '''));']); end
    catch ME
        warning('Parallel setup failed: %s. Falling back to serial.',ME.message);
        use_parallel = false;
    end
end

fprintf('\n=================================================================\n');
fprintf(' EXPERIMENT 2: OPTIMIZED UNIFIED RET-RIS ISAC\n');
fprintf(' Mode=%s, MC=%d, M=%d, N=%d, K=%d, T=%d\n', ...
    upper(RUN_MODE),ITER,M,N,K,T);
fprintf(' gamma=%.1f dB, P_total=%g:%g:%g dBm (%.3f--%.3f W)\n', ...
    gamma_dB,range_P_dBm(1),range_P_dBm(2)-range_P_dBm(1),range_P_dBm(end), ...
    range_P(1),range_P(end));
fprintf(' Road targets: '); fprintf('(%.0f,%.0f) ',[cfg.target_x;cfg.target_y]); fprintf('m\n');
fprintf(' Same M/N/K/T, road, fading, RET and seeds as Experiment 1.\n');
fprintf('=================================================================\n\n');

%% ======================== Result Buffers =================================
S=5;
Obj = nan(ITER,numP,S);
Pcomm = nan(ITER,numP,S);
Pradar = nan(ITER,numP,S);
Feasible = false(ITER,numP,S);
Theta = nan(ITER,numP,S);
Obj_fixed_alpha = nan(ITER,numP,numAlpha);
Fea_fixed_alpha = false(ITER,numP,numAlpha);
Obj_radar_upper = nan(ITER,numP);

%% ======================== Monte Carlo ====================================
fprintf('Starting Monte Carlo (%d iterations, %d total-power points)...\n', ...
    ITER, numP);

% Independent timer handles prevent Manopt/CVX internal tic/toc calls from
% interfering with the experiment timing.
t_start = tic;

% Progress is measured by completed (Monte Carlo sample, total-power point)
% pairs.  In PAPER mode this gives ITER*numP progress steps.
progressQueue = [];
if SHOW_PROGRESS
    try
        total_progress_steps = ITER * numP;
        progressQueue = parallel.pool.DataQueue;
        afterEach(progressQueue, @exp2_progress_bar);
        exp2_progress_bar(struct('action','reset','total',total_progress_steps));
    catch ME
        warning('Progress bar initialization failed: %s', ME.message);
        progressQueue = [];
    end
end

if use_parallel
    parfor iter = 1:ITER
        rng(cfg.seed_base + iter*cfg.seed_stride,'twister');
        seed=cfg.seed_base + iter*cfg.seed_stride;
        out = run_one_mc_isac(cfg,range_P,gamma,seed,func_path,progressQueue);
        Obj(iter,:,:) = reshape(out.Obj,[1,numP,S]);
        Pcomm(iter,:,:) = reshape(out.Pcomm,[1,numP,S]);
        Pradar(iter,:,:) = reshape(out.Pradar,[1,numP,S]);
        Feasible(iter,:,:) = reshape(out.Feasible,[1,numP,S]);
        Theta(iter,:,:) = reshape(out.Theta,[1,numP,S]);
        Obj_fixed_alpha(iter,:,:) = reshape(out.Obj_fixed_alpha,[1,numP,numAlpha]);
        Fea_fixed_alpha(iter,:,:) = reshape(out.Fea_fixed_alpha,[1,numP,numAlpha]);
        Obj_radar_upper(iter,:) = out.Obj_radar_upper;
    end
else
    for iter = 1:ITER
        rng(cfg.seed_base + iter*cfg.seed_stride,'twister');
        seed=cfg.seed_base + iter*cfg.seed_stride;
        out = run_one_mc_isac(cfg,range_P,gamma,seed,func_path,progressQueue);
        Obj(iter,:,:) = reshape(out.Obj,[1,numP,S]);
        Pcomm(iter,:,:) = reshape(out.Pcomm,[1,numP,S]);
        Pradar(iter,:,:) = reshape(out.Pradar,[1,numP,S]);
        Feasible(iter,:,:) = reshape(out.Feasible,[1,numP,S]);
        Theta(iter,:,:) = reshape(out.Theta,[1,numP,S]);
        Obj_fixed_alpha(iter,:,:) = reshape(out.Obj_fixed_alpha,[1,numP,numAlpha]);
        Fea_fixed_alpha(iter,:,:) = reshape(out.Fea_fixed_alpha,[1,numP,numAlpha]);
        Obj_radar_upper(iter,:) = out.Obj_radar_upper;
    end
end
runtime = toc(t_start);

%% ======================== Statistics =====================================
% Main five-scheme curves use exactly the same realizations at each power.
CommonMain=all(Feasible,3);
CommonMainRate=100*mean(CommonMain,1).';
Fea_rate=squeeze(mean(Feasible,1));

[Obj_dB,Obj_lo_dB,Obj_hi_dB,Obj_n,Obj_avg]= ...
    masked_mean_ci_to_db(Obj,CommonMain,1);
[Pcomm_dBm,Pcomm_lo_dBm,Pcomm_hi_dBm,Pcomm_n,Pcomm_avg]= ...
    masked_mean_ci_to_db(Pcomm,CommonMain,1e-3);
[Pradar_avg,Pradar_lo,Pradar_hi,Pradar_n]= ...
    masked_mean_ci_linear(Pradar,CommonMain);
[Theta_avg,Theta_lo,Theta_hi,Theta_n]= ...
    masked_mean_ci_linear(Theta,CommonMain);

upper_mask=CommonMain & isfinite(Obj_radar_upper);
[Obj_upper_dB,Obj_upper_lo_dB,Obj_upper_hi_dB,Obj_upper_n]= ...
    masked_mean_ci_to_db(Obj_radar_upper,upper_mask,1);
DeltaPradar_mW=1e3*(Pradar_avg-Pradar_avg(:,1));
RadarGain_dB=nan(numP,S); RadarGain_lo=nan(numP,S); RadarGain_hi=nan(numP,S);
RadarGain_n=zeros(numP,S);
RadarGain_dB(:,1)=0; RadarGain_lo(:,1)=0; RadarGain_hi(:,1)=0;
for ss=2:S
    [RadarGain_dB(:,ss),RadarGain_lo(:,ss),RadarGain_hi(:,ss),RadarGain_n(:,ss)]= ...
        paired_db_statistics(Obj(:,:,ss),Obj(:,:,1), ...
        Feasible(:,:,ss)&Feasible(:,:,1));
end

% Fixed-power-split comparison: all displayed curves share one common mask.
[~,a50]=min(abs(cfg.alpha_fixed_list-0.5));
[~,a70]=min(abs(cfg.alpha_fixed_list-0.7));
CommonAllocation=Feasible(:,:,1)&Feasible(:,:,5)& ...
    Fea_fixed_alpha(:,:,a50)&Fea_fixed_alpha(:,:,a70);
CommonAllocationRate=100*mean(CommonAllocation,1).';
ObjAllocation=cat(3,Obj(:,:,1),Obj(:,:,5), ...
    Obj_fixed_alpha(:,:,a50),Obj_fixed_alpha(:,:,a70));
[ObjAllocation_dB,ObjAllocation_lo,ObjAllocation_hi,ObjAllocation_n]= ...
    masked_mean_ci_to_db(ObjAllocation,CommonAllocation,1);
Obj_fixed_dB=nan(numP,numAlpha); Obj_fixed_lo=Obj_fixed_dB; Obj_fixed_hi=Obj_fixed_dB;
Obj_fixed_n=zeros(numP,numAlpha);
FixedFea_rate=squeeze(100*mean(Fea_fixed_alpha,1));
for aa=1:numAlpha
    mask_a=Feasible(:,:,1)&Feasible(:,:,5)&Fea_fixed_alpha(:,:,aa);
    [Obj_fixed_dB(:,aa),Obj_fixed_lo(:,aa),Obj_fixed_hi(:,aa),Obj_fixed_n(:,aa)]= ...
        masked_mean_ci_to_db(Obj_fixed_alpha(:,:,aa),mask_a,1);
end

[radar_targets_dB,Pmin_avg_W,Pmin_fea_rate,Pmin_samples_W, ...
    Pmin_common_rate,Pmin_common_mask,Pmin_common_n]= ...
    compute_min_power_for_snr(Obj,range_P,cfg.radar_snr_target_dB_exp2, ...
    cfg.num_radar_snr_targets);
Pmin_avg_dBm=watt2dbm(Pmin_avg_W);
labels={'Fixed RET, No RIS','Fixed RET, RIS-MinPc', ...
        'RET-NoRIS-MinPc','RET-RIS-MinPc','RET-RIS-RadarSNR'};
markers={'-o','-^','-d','-s','-p'};

if SAVE_MAIN_FIGURE
%% ======================== Main Paper Composite ===========================
% One 2x2 figure carries the main Experiment-2 evidence. The remaining
% single-panel figures are saved as supplementary diagnostics.
fig00 = figure('Position',[30 30 1450 900], 'Color', 'w');
tl00 = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

nexttile;
for s=1:S
    errorbar(range_P_dBm,Obj_dB(:,s),Obj_dB(:,s)-Obj_lo_dB(:,s), ...
        Obj_hi_dB(:,s)-Obj_dB(:,s),markers{s},'LineWidth',1.5,'MarkerSize',7); hold on;
end
plot(range_P_dBm,Obj_upper_dB,':p','LineWidth',1.4,'MarkerSize',8);
grid on; xlabel('Total transmit power P_T (dBm)'); ylabel('Weighted radar SNR (dB)');
title(sprintf('(a) Common-feasible sensing performance (min n=%d)',min(Obj_n(:))));
legend([labels,{'Radar-only reference (same RET-RIS)'}],'Location','best');

nexttile;
for s=1:S
    errorbar(range_P_dBm,Pcomm_dBm(:,s),Pcomm_dBm(:,s)-Pcomm_lo_dBm(:,s), ...
        Pcomm_hi_dBm(:,s)-Pcomm_dBm(:,s),markers{s},'LineWidth',1.5,'MarkerSize',7); hold on;
end
grid on; xlabel('Total transmit power P_T (dBm)'); ylabel('Communication power (dBm)');
title('(b) Common-feasible communication cost at fixed QoS');

nexttile;
plot(range_P_dBm,ObjAllocation_dB(:,1),'--o','LineWidth',1.4,'MarkerSize',7); hold on;
plot(range_P_dBm,ObjAllocation_dB(:,2),'-p','LineWidth',1.8,'MarkerSize',8);
plot(range_P_dBm,ObjAllocation_dB(:,3),'-.^','LineWidth',1.4,'MarkerSize',7);
plot(range_P_dBm,ObjAllocation_dB(:,4),':v','LineWidth',1.5,'MarkerSize',7);
grid on; xlabel('Total transmit power P_T (dBm)'); ylabel('Weighted radar SNR (dB)');
title(sprintf('(c) Full-budget fixed splits; common rate %.0f--%.0f%%', ...
    min(CommonAllocationRate),max(CommonAllocationRate)));
legend({'Fixed RET/no RIS','RET-RIS-RadarSNR', ...
        sprintf('Fixed radar fraction %.1f',cfg.alpha_fixed_list(a50)), ...
        sprintf('Fixed radar fraction %.1f',cfg.alpha_fixed_list(a70))}, ...
        'Location','best');

nexttile;
for s=1:S
    plot(radar_targets_dB,Pmin_avg_dBm(:,s),markers{s},'LineWidth',1.5,'MarkerSize',7); hold on;
end
grid on; xlabel('Required weighted radar SNR \Gamma_r (dB)');
ylabel('Minimum total power P_T^{min} (dBm)');
title(sprintf('(d) Equal-demand power on common reachable samples (min n=%d)', ...
    min(Pmin_common_n(:))));

sgtitle(tl00,sprintf('Communication-Sensing Power Coupling (QoS threshold = %.1f dB)',gamma_dB));
if SAVE_MAIN_FIGURE; save_figure_auto(fig00, figure_dir, 0, 'Main_Power_Coupling_Composite'); end

end

if SAVE_SUPP_FIGURES
%% ======================== Figure 1: Radar SNR =============================
fig01 = figure('Position',[80 80 900 600], 'Color', 'w');
for s=1:S
    plot(range_P_dBm,Obj_dB(:,s),markers{s},'LineWidth',1.6,'MarkerSize',8); hold on;
end
plot(range_P_dBm,Obj_upper_dB,':p','LineWidth',1.5,'MarkerSize',9);
xlabel('Total transmit power P_T (dBm)'); ylabel('Weighted radar SNR (dB)');
legend([labels,{'Radar-only reference (same RET-RIS config)'}],'Location','northwest');
grid on; title('RET-RIS ISAC Radar Performance');
if SAVE_SUPP_FIGURES; save_figure_auto(fig01, figure_dir, 1, 'S1_Weighted_Radar_SNR'); end

%% ======================== Figure 2: Visible Reallocation Gain =============
fig02 = figure('Position',[100 100 900 600], 'Color', 'w');
for s=2:S
    plot(range_P_dBm,DeltaPradar_mW(:,s),markers{s},'LineWidth',1.6,'MarkerSize',8); hold on;
end
xlabel('Total transmit power P_T (dBm)');
ylabel('\Delta P_{radar} over fixed/no-RIS baseline (mW)');
legend(labels(2:S),'Location','best'); grid on;
title('Change in Residual Sensing Power Relative to Fixed-RET/No-RIS');
yline(0,'--','HandleVisibility','off');
if SAVE_SUPP_FIGURES; save_figure_auto(fig02, figure_dir, 2, 'S2_Residual_Radar_Power_Change'); end

%% ======================== Figure 3: Communication Power ===================
fig03 = figure('Position',[120 120 900 600], 'Color', 'w');
for s=1:S
    plot(range_P_dBm,Pcomm_dBm(:,s),markers{s},'LineWidth',1.6,'MarkerSize',8); hold on;
end
xlabel('Total transmit power P_T (dBm)'); ylabel('Minimum communication power P_c (dBm)');
legend(labels,'Location','best'); grid on;
title(sprintf('Communication QoS Fixed at %.1f dB',gamma_dB));
if SAVE_SUPP_FIGURES; save_figure_auto(fig03, figure_dir, 3, 'S3_Minimum_Communication_Power'); end

%% ======================== Figure 4: Mechanism Ablation ====================
fig04 = figure('Position',[140 140 900 600], 'Color', 'w');
plot(range_P_dBm,ObjAllocation_dB(:,1),'--o','LineWidth',1.5,'MarkerSize',7); hold on;
plot(range_P_dBm,ObjAllocation_dB(:,2),'-p','LineWidth',1.9,'MarkerSize',8);
plot(range_P_dBm,ObjAllocation_dB(:,3),'-.^','LineWidth',1.5,'MarkerSize',7);
plot(range_P_dBm,ObjAllocation_dB(:,4),':v','LineWidth',1.6,'MarkerSize',7);
legend({'Dynamic: fixed tilt/no RIS','Dynamic: RET-RIS-RadarSNR', ...
        sprintf('RET+RIS fixed radar fraction %.1f',cfg.alpha_fixed_list(a50)), ...
        sprintf('RET+RIS fixed radar fraction %.1f',cfg.alpha_fixed_list(a70))}, ...
        'Location','northwest');
xlabel('Total transmit power P_T (dBm)'); ylabel('Weighted radar SNR (dB)');
grid on; title('Full-Budget Fixed Split versus Dynamic Allocation (Common Samples)');
if SAVE_SUPP_FIGURES; save_figure_auto(fig04, figure_dir, 4, 'S4_Dynamic_and_Fixed_Power_Allocation'); end

%% ======================== Figure 5: RET Tilt ===============================
fig05 = figure('Position',[160 160 1200 430], 'Color', 'w');
subplot(1,3,1); plot(range_P_dBm,Theta_avg(:,3),'-d','LineWidth',1.5); grid on;
xlabel('P_T (dBm)'); ylabel('Optimal RET tilt (deg)'); title('RET-NoRIS-MinPc');
subplot(1,3,2); plot(range_P_dBm,Theta_avg(:,4),'-s','LineWidth',1.5); grid on;
xlabel('P_T (dBm)'); ylabel('Optimal RET tilt (deg)'); title('RET-RIS-MinPc');
subplot(1,3,3); plot(range_P_dBm,Theta_avg(:,5),'-p','LineWidth',1.5); grid on;
xlabel('P_T (dBm)'); ylabel('Optimal RET tilt (deg)'); title('RET-RIS-RadarSNR');
ax=findall(fig05,'Type','axes'); set(ax,'YLim',[0 35]);
if SAVE_SUPP_FIGURES; save_figure_auto(fig05, figure_dir, 5, 'S5_Optimal_RET_Tilt_Angles'); end

%% ======================== Figure 6: Equal Sensing, Minimum Power ===========
fig06 = figure('Position',[180 180 900 600], 'Color', 'w');
for s=1:S
    plot(radar_targets_dB,Pmin_avg_dBm(:,s),markers{s},'LineWidth',1.6,'MarkerSize',8); hold on;
end
xlabel('Required weighted radar SNR \Gamma_r (dB)');
ylabel('Minimum total transmit power P_T^{min} (dBm)');
legend(labels,'Location','northwest'); grid on;
title(sprintf('Equal Communication QoS and Equal Sensing Requirement (QoS = %.1f dB)',gamma_dB));
if SAVE_SUPP_FIGURES; save_figure_auto(fig06, figure_dir, 6, 'S6_Equal_Sensing_Minimum_Total_Power'); end

%% ======================== Figure 7: Fixed-Split Feasibility ==============
fig07=figure('Position',[200 200 940 600],'Color','w'); hold on;
for aa=1:numAlpha
    plot(range_P_dBm,FixedFea_rate(:,aa),'-o','LineWidth',1.5,'MarkerSize',7);
end
plot(range_P_dBm,CommonAllocationRate,'--k','LineWidth',1.8);
grid on; ylim([0 105]); xlabel('Total transmit power P_T (dBm)');
ylabel('Feasibility rate (%)');
fixed_split_labels=cell(1,numAlpha);
for aa=1:numAlpha
    fixed_split_labels{aa}=sprintf('Fixed radar fraction %.1f',cfg.alpha_fixed_list(aa));
end
legend([fixed_split_labels,{'Common samples in allocation comparison'}], ...
    'Location','southeast');
title('Feasibility of Full-Budget Fixed Power Splits');
if SAVE_SUPP_FIGURES; save_figure_auto(fig07,figure_dir,7,'S7_Fixed_Split_Feasibility'); end


end

%% ======================== Tables / Save ==================================
fprintf('\nAverage weighted radar SNR (dB):\n');
disp(array2table([range_P_dBm(:),Obj_dB,Obj_upper_dB(:)],'VariableNames', ...
    {'P_total_dBm','Fixed_NoRIS','Fixed_RIS_MinPc','RET_NoRIS_MinPc','RET_RIS_MinPc','RET_RIS_RadarSNR','RadarUpper'}));
fprintf('\nAdditional radar power relative to fixed/no-RIS (mW):\n');
disp(array2table([range_P_dBm(:),DeltaPradar_mW(:,2:S)],'VariableNames', ...
    {'P_total_dBm','Fixed_RIS_MinPc','RET_NoRIS_MinPc','RET_RIS_MinPc','RET_RIS_RadarSNR'}));
fprintf('\nRadar-SNR gain relative to fixed/no-RIS (dB):\n');
disp(array2table([range_P_dBm(:),RadarGain_dB(:,2:S)],'VariableNames', ...
    {'P_total_dBm','Fixed_RIS_MinPc','RET_NoRIS_MinPc','RET_RIS_MinPc','RET_RIS_RadarSNR'}));
fprintf('\nMinimum total power at equal radar SNR (dBm):\n');
disp(array2table([radar_targets_dB(:),Pmin_avg_dBm],'VariableNames', ...
    {'RadarSNR_Target','Fixed_NoRIS','Fixed_RIS_MinPc','RET_NoRIS_MinPc','RET_RIS_MinPc','RET_RIS_RadarSNR'}));
fprintf('\nEqual-SNR minimum-power feasibility rate:\n');
disp(array2table([radar_targets_dB(:),Pmin_fea_rate],'VariableNames', ...
    {'RadarSNR_Target','Fixed_NoRIS','Fixed_RIS_MinPc','RET_NoRIS_MinPc','RET_RIS_MinPc','RET_RIS_RadarSNR'}));


fprintf('\nCommon-feasible rate for the five main schemes (%%):\n');
disp(array2table([range_P_dBm(:),CommonMainRate], ...
    'VariableNames',{'P_total_dBm','CommonMain_percent'}));
fprintf('\nFixed-split feasibility and common comparison rate (%%):\n');
disp(array2table([range_P_dBm(:),FixedFea_rate,CommonAllocationRate], ...
    'VariableNames',{'P_total_dBm','Alpha_03','Alpha_05','Alpha_07','CommonAllocation_percent'}));
fprintf('\nEqual-SNR common reachability rate (%%):\n');
disp(array2table([radar_targets_dB(:),Pmin_common_rate], ...
    'VariableNames',{'RadarSNR_Target','CommonReachable_percent'}));

results_dir=fullfile(root,'results'); if ~exist(results_dir,'dir'); mkdir(results_dir); end
save(fullfile(results_dir,sprintf('exp2_revised_%s.mat',lower(RUN_MODE))), ...
    'cfg','range_P','range_P_dBm','Obj','Pcomm','Pradar','Feasible','Theta', ...
    'Obj_fixed_alpha','Fea_fixed_alpha','Obj_radar_upper','CommonMain', ...
    'CommonMainRate','Obj_dB','Obj_lo_dB','Obj_hi_dB','Obj_n', ...
    'Pcomm_dBm','Pcomm_lo_dBm','Pcomm_hi_dBm','Pcomm_n','Pradar_avg', ...
    'CommonAllocation','CommonAllocationRate','ObjAllocation_dB', ...
    'ObjAllocation_lo','ObjAllocation_hi','ObjAllocation_n','FixedFea_rate', ...
    'RadarGain_dB','RadarGain_lo','RadarGain_hi','RadarGain_n', ...
    'Obj_upper_dB','Obj_upper_lo_dB','Obj_upper_hi_dB','Obj_upper_n', ...
    'radar_targets_dB','Pmin_avg_W','Pmin_fea_rate','Pmin_samples_W', ...
    'Pmin_common_rate','Pmin_common_mask','Pmin_common_n');
Tmain=array2table([range_P_dBm(:),CommonMainRate,Obj_dB,Pcomm_dBm], ...
    'VariableNames',{'P_total_dBm','CommonMain_percent', ...
    'Obj1_dB','Obj2_dB','Obj3_dB','Obj4_dB','Obj5_dB', ...
    'Pc1_dBm','Pc2_dBm','Pc3_dBm','Pc4_dBm','Pc5_dBm'});
Talloc=array2table([range_P_dBm(:),CommonAllocationRate,FixedFea_rate,ObjAllocation_dB], ...
    'VariableNames',{'P_total_dBm','CommonAllocation_percent', ...
    'FeaAlpha03_percent','FeaAlpha05_percent','FeaAlpha07_percent', ...
    'DynamicFixedNoRIS_dB','DynamicProposed_dB','FixedAlpha05_dB','FixedAlpha07_dB'});
writetable(Tmain,fullfile(results_dir,sprintf('exp2_main_statistics_%s.csv',lower(RUN_MODE))));
writetable(Talloc,fullfile(results_dir,sprintf('exp2_allocation_statistics_%s.csv',lower(RUN_MODE))));
total_time = toc(t_total);
fprintf('\n');
fprintf('==================================================\n');
fprintf('  SIMULATION COMPLETED (Exp2 - RET-RIS ISAC)\n');
fprintf('  ------------------------------------------------\n');
fprintf('  Monte Carlo time:  %8.1f s  (%6.2f min / %5.2f h)\n', ...
    runtime, runtime/60, runtime/3600);
fprintf('  Total runtime:     %8.1f s  (%6.2f min / %5.2f h)\n', ...
    total_time, total_time/60, total_time/3600);
fprintf('  Avg per MC run:    %8.1f s  (%6.2f min)\n', ...
    runtime/ITER, runtime/ITER/60);
fprintf('  Avg per power task:%8.1f s\n', runtime/(ITER*numP));
fprintf('  Mode:              %s\n', upper(RUN_MODE));
if SAVE_FIGURES; fprintf('  Figure folder:     %s\n',figure_dir); else; fprintf('  Per-exp figures:   not exported (%s profile)\n',FIG_POLICY.profile); end
if use_parallel
    fprintf('  Workers:           %d\n', num_workers);
end
fprintf('==================================================\n');

%% ========================================================================
%% Local Functions
%% ========================================================================

function out = run_one_mc_isac(cfg,range_P,gamma,seed,func_path,progressQueue)
% One paired roadside realization.  The five schemes are evaluated only
% through the shared evaluator used by Exp. 1 and Exp. 3.
if ~isempty(func_path) && exist(func_path,'dir'); addpath(func_path); end
scen=generate_isac_scenario(cfg,cfg.N,cfg.K,seed,[],[],struct());
solver_seed=seed+cfg.solver_seed_offset_sweep;
sw=evaluate_ret_ris_sweep(cfg,scen,gamma,range_P,solver_seed, ...
    struct('use_continuation',cfg.exp23_use_continuation, ...
    'seed_policy','operating_point'));

numP=numel(range_P); numAlpha=numel(cfg.alpha_fixed_list);
out.Obj=sw.Obj; out.Pcomm=sw.Pcomm; out.Pradar=sw.Pradar;
out.Feasible=sw.Feasible; out.Theta=sw.Theta;
out.Obj_fixed_alpha=nan(numP,numAlpha);
out.Fea_fixed_alpha=false(numP,numAlpha);
out.Obj_radar_upper=nan(1,numP);

for pIdx=1:numP
    if sw.Feasible(pIdx,5) && ~isempty(sw.Phi{pIdx,5})
        theta5=sw.Theta(pIdx,5); phi5=sw.Phi{pIdx,5};
        for a=1:numAlpha
            rf=evaluate_fixed_split_configuration(cfg,scen,theta5,phi5, ...
                gamma,range_P(pIdx),cfg.alpha_fixed_list(a),true);
            if rf.feasible
                out.Obj_fixed_alpha(pIdx,a)=rf.obj_snr;
                out.Fea_fixed_alpha(pIdx,a)=true;
            end
        end
        out.Obj_radar_upper(pIdx)=evaluate_radar_only_reference(cfg,scen, ...
            theta5,phi5,range_P(pIdx),true);
    end
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


function [snr_targets_dB,Pmin_avg_W,Fea_rate,Pmin_samples_W, ...
    CommonRate,CommonMask,CommonN]=compute_min_power_for_snr(Obj,range_P,target_cfg,num_targets)
% Convert radar-SNR-vs-power curves into the inverse experiment: minimum
% total power required to meet the same radar SNR target.
[ITER,~,S]=size(Obj);
range_P=range_P(:).';
Obj_dB_samples=nan(size(Obj));
valid=isfinite(Obj)&Obj>0;
Obj_dB_samples(valid)=10*log10(Obj(valid));

if nargin<4 || isempty(num_targets); num_targets=5; end
if isempty(target_cfg)
    avg_dB=nan(numel(range_P),S);
    common_power=all(isfinite(Obj)&Obj>0,3);
    for ss=1:S
        [tmp_db,~,~,~]=masked_mean_ci_to_db(Obj(:,:,ss),common_power,1);
        avg_dB(:,ss)=tmp_db;
    end
    lo_s=nan(1,S); hi_s=nan(1,S);
    for s=1:S
        v=avg_dB(:,s); v=v(isfinite(v));
        if ~isempty(v); lo_s(s)=min(v); hi_s(s)=max(v); end
    end
    lo=max(lo_s(isfinite(lo_s))); hi=min(hi_s(isfinite(hi_s)));
    if isempty(lo) || isempty(hi) || hi<=lo
        v=avg_dB(:,1); v=v(isfinite(v));
        lo=min(v); hi=max(v);
    end
    span=max(hi-lo,0.2); margin=min(0.08*span,0.2);
    lo=lo+margin; hi=hi-margin;
    if hi<=lo; hi=lo+0.2; end
    snr_targets_dB=linspace(lo,hi,num_targets);
else
    snr_targets_dB=target_cfg(:).';
end

Pmin_samples_W=nan(ITER,numel(snr_targets_dB),S);
for it=1:ITER
    for s=1:S
        snr=squeeze(Obj_dB_samples(it,:,s));
        snr=snr(:).';
        good=isfinite(snr)&isfinite(range_P);
        p=range_P(good); snr=snr(good);
        if isempty(snr); continue; end
        [p,ord]=sort(p); snr=snr(ord);
        % Do not enforce a monotone envelope. The first sampled power that
        % reaches the target is used, with interpolation only across its
        % immediately preceding below-target sample. Any genuine numerical
        % non-monotonicity remains visible in the saved raw curves.
        for q=1:numel(snr_targets_dB)
            target=snr_targets_dB(q);
            idx=find(snr>=target,1,'first');
            if isempty(idx); continue; end
            if idx==1 || snr(idx)<=snr(idx-1)+1e-12
                p_req=p(idx);
            else
                frac=(target-snr(idx-1))/(snr(idx)-snr(idx-1));
                frac=min(max(frac,0),1);
                p_req=p(idx-1)+frac*(p(idx)-p(idx-1));
            end
            Pmin_samples_W(it,q,s)=p_req;
        end
    end
end

Q=numel(snr_targets_dB);
Pmin_avg_W=nan(Q,S);
Fea_rate=nan(Q,S);
for ss=1:S
    for qq=1:Q
        Fea_rate(qq,ss)=mean(isfinite(Pmin_samples_W(:,qq,ss)));
    end
end
CommonMask=all(isfinite(Pmin_samples_W),3);
CommonRate=100*mean(CommonMask,1).';
CommonN=sum(CommonMask,1).';
for qq=1:Q
    ok=CommonMask(:,qq);
    if ~any(ok); continue; end
    for ss=1:S
        Pmin_avg_W(qq,ss)=mean(Pmin_samples_W(ok,qq,ss));
    end
end
end


%% ======================== Figure Saving Helper ============================
function save_figure_auto(fig_handle, figure_dir, figure_index, figure_title)
% SAVE_FIGURE_AUTO Save one figure with a sequential, descriptive filename.
% Each figure is saved as:
%   01_Name.png  - 300 dpi image for reports and papers
%   01_Name.fig  - editable MATLAB figure source

    if ~exist(figure_dir, 'dir')
        mkdir(figure_dir);
    end

    safe_title = regexprep(figure_title, '[^A-Za-z0-9_\-]', '_');
    base_name = sprintf('%02d_%s', figure_index, safe_title);
    png_path = fullfile(figure_dir, [base_name, '.png']);
    fig_path = fullfile(figure_dir, [base_name, '.fig']);

    drawnow;
    try
        exportgraphics(fig_handle, png_path, 'Resolution', 300);
    catch
        % Compatibility fallback for MATLAB releases without exportgraphics.
        print(fig_handle, png_path, '-dpng', '-r300');
    end

    try
        savefig(fig_handle, fig_path);
    catch ME
        warning('Could not save editable FIG file %s: %s', fig_path, ME.message);
    end

    fprintf('Saved Figure %02d: %s\n', figure_index, png_path);
end

function exp2_progress_bar(message)
% EXP2_PROGRESS_BAR
% Single-line command-window progress bar with elapsed time and ETA.
% The callback runs on the MATLAB client even when workers call send().

    persistent total_steps completed_steps start_timer;

    % Initialization/reset message.
    if isstruct(message)
        if isfield(message,'action') && strcmp(message.action,'reset')
            total_steps = double(message.total);
            completed_steps = 0;
            start_timer = tic;
            fprintf(['Progress: [                              ] ', ...
                '0/%d (  0.0%%)  Elapsed 00:00:00  ETA --:--:--'], ...
                total_steps);
        end
        return;
    end

    if isempty(total_steps) || total_steps <= 0
        return;
    end

    completed_steps = min(total_steps, completed_steps + double(message));
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
    progress_text = [repmat('=',1,filled_length), ...
        repmat(' ',1,bar_length-filled_length)];

    elapsed_text = format_exp2_time(elapsed_sec);
    if isfinite(eta_sec)
        eta_text = format_exp2_time(eta_sec);
    else
        eta_text = '--:--:--';
    end

    fprintf(['\rProgress: [%s] %d/%d (%5.1f%%)  ', ...
        'Elapsed %s  ETA %s'], ...
        progress_text, completed_steps, total_steps, percentage, ...
        elapsed_text, eta_text);
    drawnow limitrate;

    if completed_steps >= total_steps
        fprintf('\n');
        total_steps = [];
        completed_steps = [];
        start_timer = [];
    end
end

function time_text = format_exp2_time(time_seconds)
% Convert seconds to HH:MM:SS.
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

