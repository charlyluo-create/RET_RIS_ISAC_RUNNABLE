% =========================================================================
% Experiment 7 SUBMISSION-READY: Same-Model Algorithm Comparison
%
% Fair baselines under the same downlink vehicular-ISAC model:
%   1) Fixed RET, No RIS
%   2) Fixed RET, RIS-RadarSNR
%   3) Sequential RET then RIS
%   4) Random RET-RIS search (supplementary sanity baseline)
%   5) Fine AO block-coordinate phase search from the same Stage-I seed (high-compute reference)
%   6) RET-RIS-MinPc (different-objective mechanism baseline)
%   7) Accelerated Stage-I RCG-AO (pre-refinement snapshot; zero extra solve)
%   8) Proposed MA-RCG-MF v2 (sequential MF block refinement + ATR)
%
% The main paper figure uses schemes 1,2,3,5,7,8. Scheme 7 is the exact
% pre-refinement endpoint produced inside the same MA-RCG-MF execution, so
% Stage-II gain and cost are measured without rerunning Stage I. The figure
% reports common-feasible objective statistics, complete executable runtime,
% true-objective improvement and paired SNR gaps. Random search and MinPc
% remain in the exported table/supplement. The strong BCS baseline starts from
% the identical Stage-I point and is not weakened to favor the proposed method.
% =========================================================================

clear; close all; clc;
t_total=tic;
SHOW_PROGRESS=true;

%% Configuration
env_mode=getenv('RET_RIS_RUN_MODE');
if isempty(env_mode); RUN_MODE='paper'; else; RUN_MODE=lower(env_mode); end
cfg=unified_ret_ris_isac_config(RUN_MODE);
FIG_POLICY=get_figure_policy(RUN_MODE);
SAVE_MAIN_FIGURE=FIG_POLICY.save_main;
SAVE_SUPP_FIGURES=FIG_POLICY.save_supplementary;
SAVE_FIGURES=FIG_POLICY.save_any;
ITER=cfg.ITER_EXP7;
gamma_dB=cfg.gamma_exp7_dB;
gamma=10^(gamma_dB/10);
P_total_dBm=cfg.P_total_exp7_dBm;
P_total=cfg.P_total_exp7;
SolverMode=getenv('RET_RIS_COMM_SOLVER');
if isempty(SolverMode); SolverMode=cfg.comm_solver_default; setenv('RET_RIS_COMM_SOLVER',SolverMode); end
SolverMode=lower(strtrim(SolverMode));
RuntimeAccountingPolicy=cfg.comm_solver_runtime_accounting;
assert(any(strcmp(SolverMode,{'hybrid','duality','fast','cvx','socp'})), ...
    'Unsupported RET_RIS_COMM_SOLVER mode: %s',SolverMode);

labels={'Fixed RET, No RIS', ...
        'Fixed RET, RIS-RadarSNR', ...
        'Sequential RET then RIS', ...
        'Random RET-RIS search', ...
        'Fine AO block-coordinate search', ...
        'RET-RIS-MinPc', ...
        'Accelerated Stage-I RCG-AO', ...
        'Proposed MA-RCG-MF'};
S=numel(labels);
main_idx=[1 2 3 5 7 8];

root=fileparts(mfilename('fullpath'));
if isempty(root); root=pwd; end
addpath(root); addpath(fullfile(root,'function'));
required={'evaluate_ret_ris_schemes','evaluate_fixed_ret_ris_config', ...
    'optimize_fixed_ret_ris_radar','optimize_ret_ris_block_coordinate'};
for ii=1:numel(required)
    if exist(required{ii},'file')~=2; error('Missing required function: %s.m',required{ii}); end
end
if exist('cvx_begin','file')~=2
    error('CVX is required. Install CVX and run cvx_setup first.');
end
run_stamp=datestr(now,'yyyymmdd_HHMMSS');
figure_dir=fullfile(root,sprintf('Exp7_Figures_%s_%s',lower(RUN_MODE),run_stamp));
if SAVE_FIGURES && ~exist(figure_dir,'dir'); mkdir(figure_dir); end

fprintf('\n===============================================================\n');
fprintf(' Experiment 7: Same-Model Algorithm Comparison\n');
fprintf(' Mode=%s, MC=%d, gamma=%.1f dB, P_T=%.1f dBm\n', ...
    upper(RUN_MODE),ITER,gamma_dB,P_total_dBm);
fprintf(' BCS baseline: %d blocks, %d phase levels, %d sweeps\n', ...
    cfg.exp7_bcs_blocks,cfg.exp7_bcs_phase_levels,cfg.exp7_bcs_sweeps);
if isfield(cfg,'proposed_stage2_variant')
    fprintf(' Proposed Stage-II: %s, %d blocks, up to %d sequential sweeps\n', ...
        cfg.proposed_stage2_variant,cfg.v2_blocks,cfg.v2_max_sweeps);
end
fprintf(' Common communication solver backend policy: %s (identical for all schemes)\n',upper(SolverMode));
fprintf(' Runtime accounting: %s\n',RuntimeAccountingPolicy);
fprintf('===============================================================\n');

Obj=nan(ITER,S); Pcomm=nan(ITER,S); Feasible=false(ITER,S); Runtime=nan(ITER,S);
% V2 can accept one exact update per block in each sequential sweep, so the
% history buffer must accommodate many more accepted physical-objective moves
% than the legacy two-round global polish.
if isfield(cfg,'v2_blocks') && isfield(cfg,'v2_max_sweeps')
    maxHist=cfg.max_ao_iter+8+cfg.v2_blocks*cfg.v2_max_sweeps+cfg.v2_max_sweeps;
else
    Rpol=numel(cfg.proposed_polish_blocks_schedule);
    maxHist=cfg.max_ao_iter+4+2*Rpol;
end
HistRadar=nan(ITER,maxHist);
HistMinPc=nan(ITER,maxHist);
RefineGain_dB=nan(ITER,1); RefineExactEvals=nan(ITER,1); RefineAccepts=nan(ITER,1);
Stage1ExactEvals=nan(ITER,1); Stage1Accepts=nan(ITER,1); TotalExactEvals=nan(ITER,1); ProxyScreened=nan(ITER,1);
% Runtime decomposition used to audit the performance-complexity story.
Stage1Runtime=nan(ITER,1); Stage2Runtime=nan(ITER,1);
BCSRefineRuntime=nan(ITER,1); BCSExactEvals=nan(ITER,1);
% V2 Stage-II audit columns:
% [sweeps, block_visits, budget_expansions, ret_exact, ris_exact,
%  mean_final_trust_deg, max_final_trust_deg, accepted_block_updates]
V2Metrics=nan(ITER,8);
SolverTelemetry=cell(ITER,1); % each cell is S-by-8, see comm_solver_stats_summary

%% Parallel execution
use_parallel=true; progressQueue=[];
if SHOW_PROGRESS
    try
        progressQueue=parallel.pool.DataQueue;
        afterEach(progressQueue,@(~)exp7_progress(ITER));
    catch
        progressQueue=[];
    end
end
if use_parallel
    try
        poolobj=gcp('nocreate');
        if isempty(poolobj); poolobj=parpool(min(cfg.num_workers,ITER)); end
        pctRunOnAll(['addpath(''' root ''');']);
        pctRunOnAll(['addpath(''' fullfile(root,'function') ''');']);
        pctRunOnAll(sprintf('setenv(''RET_RIS_COMM_SOLVER'',''%s'');',SolverMode));
    catch ME
        warning('Parallel setup failed: %s. Running serially.',ME.message);
        use_parallel=false;
    end
end

runtime=tic;
if use_parallel
    parfor mc=1:ITER
        [Obj(mc,:),Pcomm(mc,:),Feasible(mc,:),Runtime(mc,:), ...
            HistRadar(mc,:),HistMinPc(mc,:),RefineGain_dB(mc), ...
            RefineExactEvals(mc),RefineAccepts(mc),Stage1ExactEvals(mc),Stage1Accepts(mc), ...
            TotalExactEvals(mc),ProxyScreened(mc),Stage1Runtime(mc),Stage2Runtime(mc), ...
            BCSRefineRuntime(mc),BCSExactEvals(mc),V2Metrics(mc,:),SolverTelemetry{mc}]= ...
            run_one_exp7(cfg,mc,gamma,P_total,maxHist);
        if ~isempty(progressQueue); send(progressQueue,1); end
    end
else
    for mc=1:ITER
        [Obj(mc,:),Pcomm(mc,:),Feasible(mc,:),Runtime(mc,:), ...
            HistRadar(mc,:),HistMinPc(mc,:),RefineGain_dB(mc), ...
            RefineExactEvals(mc),RefineAccepts(mc),Stage1ExactEvals(mc),Stage1Accepts(mc), ...
            TotalExactEvals(mc),ProxyScreened(mc),Stage1Runtime(mc),Stage2Runtime(mc), ...
            BCSRefineRuntime(mc),BCSExactEvals(mc),V2Metrics(mc,:),SolverTelemetry{mc}]= ...
            run_one_exp7(cfg,mc,gamma,P_total,maxHist);
        if SHOW_PROGRESS; exp7_progress(ITER); end
    end
end
mc_runtime=toc(runtime);

% Unpack worker-local solver telemetry. Columns are documented by
% comm_solver_stats_summary: calls, duality success, CVX attempts, hybrid
% fallbacks, fallback rate, mean FPI iterations, mean primal-dual gap,
% communication-solver wall time.
SolverInnerCalls=nan(ITER,S); SolverDualitySuccess=nan(ITER,S);
SolverCVXAttempts=nan(ITER,S); SolverFallbacks=nan(ITER,S);
SolverFallbackRate=nan(ITER,S); SolverMeanFPIIterations=nan(ITER,S);
SolverMeanPrimalDualGap=nan(ITER,S); SolverWallTime=nan(ITER,S);
for mc=1:ITER
    Z=SolverTelemetry{mc};
    if isempty(Z)||size(Z,1)~=S||size(Z,2)~=8; continue; end
    SolverInnerCalls(mc,:)=Z(:,1).'; SolverDualitySuccess(mc,:)=Z(:,2).';
    SolverCVXAttempts(mc,:)=Z(:,3).'; SolverFallbacks(mc,:)=Z(:,4).';
    SolverFallbackRate(mc,:)=Z(:,5).'; SolverMeanFPIIterations(mc,:)=Z(:,6).';
    SolverMeanPrimalDualGap(mc,:)=Z(:,7).'; SolverWallTime(mc,:)=Z(:,8).';
end

%% Statistics
% Per-scheme summaries are retained for the full exported table.
[ObjMean_dB,ObjLo_dB,ObjHi_dB,Nfea]=linear_mean_ci_to_db(Obj,Feasible);
[PcommMean_dBm,PcommLo_dBm,PcommHi_dBm,~]=linear_mean_ci_to_db(Pcomm,Feasible,1e3);
FeaRate=100*mean(Feasible,1);
RuntimeMedian=nanmedian_local(Runtime,1);
RuntimeQ1=nanquantile_local(Runtime,0.25,1);
RuntimeQ3=nanquantile_local(Runtime,0.75,1);
InnerSolveMedian=nanmedian_local(SolverInnerCalls,1);
CVXAttemptsMedian=nanmedian_local(SolverCVXAttempts,1);
FallbacksMedian=nanmedian_local(SolverFallbacks,1);
FallbackRateMedian=100*nanmedian_local(SolverFallbackRate,1);
FPIIterationsMedian=nanmedian_local(SolverMeanFPIIterations,1);
PrimalDualGapMedian=nanmedian_local(SolverMeanPrimalDualGap,1);
SolverWallMedian=nanmedian_local(SolverWallTime,1);

% The main comparison uses exactly the same channel realizations for every
% sensing-oriented algorithm, avoiding selection bias from scheme-specific
% feasible subsets.
CommonMain=all(Feasible(:,main_idx),2);
if nnz(CommonMain)<3
    warning('Only %d common-feasible samples are available for the main comparison.',nnz(CommonMain));
end
CommonMask=false(size(Feasible));
CommonMask(:,main_idx)=repmat(CommonMain,1,numel(main_idx));
[ObjCommonMean_dB,ObjCommonLo_dB,ObjCommonHi_dB,Ncommon]= ...
    linear_mean_ci_to_db(Obj,CommonMask);
[PcommCommonMean_dBm,PcommCommonLo_dBm,PcommCommonHi_dBm,~]= ...
    linear_mean_ci_to_db(Pcomm,CommonMask,1e3);

proposed_idx=8;
baseline_idx=main_idx(main_idx~=proposed_idx);
[GainMean_dB,GainLo_dB,GainHi_dB,GainN]= ...
    paired_db_gain(Obj,Feasible,proposed_idx,baseline_idx,CommonMain);

Tsummary=table(labels(:),ObjMean_dB(:),ObjLo_dB(:),ObjHi_dB(:), ...
    PcommMean_dBm(:),PcommLo_dBm(:),PcommHi_dBm(:), ...
    FeaRate(:),Nfea(:),repmat(ITER,S,1),RuntimeMedian(:),RuntimeQ1(:),RuntimeQ3(:), ...
    InnerSolveMedian(:),CVXAttemptsMedian(:),FallbacksMedian(:),FallbackRateMedian(:), ...
    FPIIterationsMedian(:),PrimalDualGapMedian(:),SolverWallMedian(:), ...
    'VariableNames',{'Scheme','RadarSNR_Mean_dB','RadarSNR_CI95_Low_dB', ...
    'RadarSNR_CI95_High_dB','Pcomm_Mean_dBm','Pcomm_CI95_Low_dBm', ...
    'Pcomm_CI95_High_dBm','Feasibility_percent','Feasible_samples', ...
    'Total_samples','Runtime_median_s','Runtime_Q1_s','Runtime_Q3_s', ...
    'InnerSolveCalls_median','CVXAttempts_median','HybridFallbacks_median', ...
    'HybridFallbackRate_percent','MeanFPIIterations_median','MeanPrimalDualGap_median', ...
    'CommSolverWallTime_median_s'});
Tpaired=table(labels(baseline_idx).',GainMean_dB(:),GainLo_dB(:),GainHi_dB(:),GainN(:), ...
    'VariableNames',{'Baseline','ProposedMinusBaseline_Mean_dB', ...
    'CI95_Low_dB','CI95_High_dB','Paired_samples'});
Trefine=table((1:ITER).',RefineGain_dB,Stage1ExactEvals,Stage1Accepts,RefineExactEvals, ...
    TotalExactEvals,ProxyScreened,RefineAccepts, ...
    'VariableNames',{'MC_index','StageII_Gain_dB','StageI_exact_evals','StageI_accepted_updates', ...
    'StageII_exact_evals','Total_exact_evals','Proxy_screened_candidates', ...
    'Accepted_refinement_updates'});
Truntime=table((1:ITER).',Stage1Runtime,Stage2Runtime,BCSRefineRuntime,BCSExactEvals, ...
    Runtime(:,5),Runtime(:,8), ...
    'VariableNames',{'MC_index','MA_complete_stageI_s','MA_stageII_refine_s', ...
    'BlockSearch_refine_only_s','BlockSearch_exact_evals','BlockSearch_complete_s', ...
    'MA_complete_s'});
Tv2=table((1:ITER).',V2Metrics(:,1),V2Metrics(:,2),V2Metrics(:,3), ...
    V2Metrics(:,4),V2Metrics(:,5),V2Metrics(:,6),V2Metrics(:,7),V2Metrics(:,8), ...
    'VariableNames',{'MC_index','Sweeps_completed','Block_visits','Budget_expansions', ...
    'RET_exact_evals','RIS_exact_evals','Mean_final_trust_deg','Max_final_trust_deg', ...
    'Accepted_block_updates'});
disp(Tsummary);
fprintf('\nCommon-feasible samples used by the main figure: %d/%d\n',nnz(CommonMain),ITER);
disp(Tpaired);
fprintf('\nMA-RCG-MF refinement gain: median %.3f dB, mean %.3f dB.\n', ...
    median(RefineGain_dB(isfinite(RefineGain_dB))),mean(RefineGain_dB(isfinite(RefineGain_dB))));
fprintf('Median exact evaluations: Stage-I %.1f, Stage-II %.1f, total %.1f; median proxy-screened candidates %.1f.\n', ...
    median(Stage1ExactEvals(isfinite(Stage1ExactEvals))), ...
    median(RefineExactEvals(isfinite(RefineExactEvals))), ...
    median(TotalExactEvals(isfinite(TotalExactEvals))), ...
    median(ProxyScreened(isfinite(ProxyScreened))));
fprintf('Proposed solver telemetry: median inner solves %.1f, CVX fallbacks %.1f, fallback rate %.3f%%.\n', ...
    InnerSolveMedian(proposed_idx),FallbacksMedian(proposed_idx),FallbackRateMedian(proposed_idx));
fprintf('BCS solver telemetry: median inner solves %.1f, CVX fallbacks %.1f, fallback rate %.3f%%.\n', ...
    InnerSolveMedian(5),FallbacksMedian(5),FallbackRateMedian(5));
fprintf('Median runtime decomposition: MA Stage-I %.3f s + Stage-II %.3f s; Block-search refinement %.3f s.\n', ...
    median(Stage1Runtime(isfinite(Stage1Runtime))),median(Stage2Runtime(isfinite(Stage2Runtime))), ...
    median(BCSRefineRuntime(isfinite(BCSRefineRuntime))));
fprintf('Median fine block-search exact evaluations in refinement: %.1f.\n', ...
    median(BCSExactEvals(isfinite(BCSExactEvals))));

%% Convergence normalization
RadarProgress=normalize_max_progress(HistRadar);
MinPcProgress=normalize_min_progress(HistMinPc);
RadarGainTrace_dB=history_gain_db(HistRadar);
[rad_gain_med,rad_gain_q1,rad_gain_q3]=column_median_iqr(RadarGainTrace_dB);
[rad_med,rad_q1,rad_q3]=column_median_iqr(RadarProgress);
[min_med,min_q1,min_q3]=column_median_iqr(MinPcProgress);
it=0:maxHist-1;
rad_updates=count_effective_updates(HistRadar,true);
min_updates=count_effective_updates(HistMinPc,false);

if SAVE_MAIN_FIGURE
%% Main paper figure: objective, runtime, convergence and paired gap
fig1=figure('Color','w','Position',[50 50 1500 900]);
tl=tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

nexttile;
x=1:numel(main_idx);
y=ObjCommonMean_dB(main_idx);
elo=y-ObjCommonLo_dB(main_idx); ehi=ObjCommonHi_dB(main_idx)-y;
errorbar(x,y,elo,ehi,'o','LineWidth',1.6,'MarkerSize',8);
grid on; xlim([0.35 numel(main_idx)+0.65]);
set(gca,'XTick',x,'XTickLabel',labels(main_idx),'XTickLabelRotation',14);
ylabel('Weighted radar SNR (dB)');
title(sprintf('(a) Common-feasible objective comparison (n = %d)',nnz(CommonMain)));

nexttile;
y=RuntimeMedian(main_idx);
elo=y-RuntimeQ1(main_idx); ehi=RuntimeQ3(main_idx)-y;
errorbar(x,y,elo,ehi,'s','LineWidth',1.6,'MarkerSize',8);
grid on; xlim([0.35 numel(main_idx)+0.65]);
set(gca,'XTick',x,'XTickLabel',labels(main_idx),'XTickLabelRotation',14);
ylabel('Median cumulative runtime (s)');
if all(y(isfinite(y))>0); set(gca,'YScale','log'); end
title('(b) Complete executable runtime');

nexttile;
plot_band(it,rad_gain_med,rad_gain_q1,rad_gain_q3,'o');
grid on;
stage1_boundary=median(Stage1Accepts(isfinite(Stage1Accepts)));
if ~isempty(stage1_boundary) && isfinite(stage1_boundary)
    xline(stage1_boundary,'--','Stage I end','LabelVerticalAlignment','bottom');
end
xlabel('Accepted update index');
ylabel('Radar-SNR gain over initialization (dB)');
title(sprintf('(c) True-objective improvement incl. refinement (median updates = %.1f)',nanmedian_scalar(rad_updates)));

nexttile;
xg=1:(numel(baseline_idx)+1); yg=[GainMean_dB(:).' 0];
lo=[GainLo_dB(:).' 0]; hi=[GainHi_dB(:).' 0];
elo=yg-lo; ehi=hi-yg;
errorbar(xg,yg,elo,ehi,'d','LineWidth',1.6,'MarkerSize',8); hold on;
yline(0,'--k','Equal objective','HandleVisibility','off');
grid on; xlim([0.5 numel(xg)+0.5]);
gapLabels=[labels(baseline_idx) {'MA-RCG-MF (reference)'}];
set(gca,'XTick',xg,'XTickLabel',gapLabels,'XTickLabelRotation',14);
ylabel('MA-RCG-MF minus comparator radar SNR (dB)');
title('(d) Paired radar-SNR gain with 95% CI');

sgtitle(tl,sprintf('Algorithm Validation (QoS = %.1f dB, P_T = %.1f dBm)', ...
    gamma_dB,P_total_dBm));
if SAVE_MAIN_FIGURE; save_exp7_figure(fig1,figure_dir,'01_Main_Algorithm_Validation'); end

end

if SAVE_SUPP_FIGURES
%% Supplementary all-scheme objective figure
fig2=figure('Color','w','Position',[100 100 1350 620]);
xall=1:S; y=ObjMean_dB; elo=y-ObjLo_dB; ehi=ObjHi_dB-y;
errorbar(xall,y,elo,ehi,'o','LineWidth',1.6,'MarkerSize',8);
grid on; xlim([0.5 S+0.5]);
set(gca,'XTick',xall,'XTickLabel',labels,'XTickLabelRotation',18);
ylabel('Weighted radar SNR (dB)');
title('Supplementary: all same-model baselines');
if SAVE_SUPP_FIGURES; save_exp7_figure(fig2,figure_dir,'S1_All_Algorithm_Baselines'); end

%% Supplementary communication-cost figure on the common main sample set
fig3=figure('Color','w','Position',[120 120 1200 620]);
y=PcommCommonMean_dBm(main_idx);
elo=y-PcommCommonLo_dBm(main_idx); ehi=PcommCommonHi_dBm(main_idx)-y;
errorbar(x,y,elo,ehi,'s','LineWidth',1.6,'MarkerSize',8);
grid on; xlim([0.35 numel(main_idx)+0.65]);
set(gca,'XTick',x,'XTickLabel',labels(main_idx),'XTickLabelRotation',14);
ylabel('Communication power (dBm)');
title(sprintf('Supplementary: communication cost on common-feasible samples (n = %d)',nnz(CommonMain)));
if SAVE_SUPP_FIGURES; save_exp7_figure(fig3,figure_dir,'S2_Common_Feasible_Communication_Cost'); end


%% Supplementary performance-runtime Pareto map (all algorithms)
fig4=figure('Color','w','Position',[140 140 1050 700]); hold on;
for ss=1:S
    if isfinite(RuntimeMedian(ss)) && isfinite(ObjMean_dB(ss))
        scatter(RuntimeMedian(ss),ObjMean_dB(ss),110,'o','LineWidth',1.7);
        text(RuntimeMedian(ss)*1.03,ObjMean_dB(ss),labels{ss}, ...
            'FontSize',9,'VerticalAlignment','middle');
    end
end
grid on; box on; xlabel('Median complete runtime (s)');
ylabel('Mean weighted radar SNR (dB)');
if all(RuntimeMedian(isfinite(RuntimeMedian))>0); set(gca,'XScale','log'); end
title('Performance-Runtime Tradeoff across All Same-Model Algorithms');
if SAVE_SUPP_FIGURES; save_exp7_figure(fig4,figure_dir,'S3_Performance_Runtime_Pareto'); end

%% Supplementary normalized convergence (diagnostic only)
fig5=figure('Color','w','Position',[160 160 900 600]);
plot_band(it,rad_med,rad_q1,rad_q3,'o'); grid on; ylim([-0.05 1.05]);
xlabel('AO outer iteration'); ylabel('Normalized objective progress');
title('Normalized Convergence Diagnostic (Not an Optimality Certificate)');
if SAVE_SUPP_FIGURES; save_exp7_figure(fig5,figure_dir,'S4_Normalized_Convergence_Diagnostic'); end

end

%% Save exact data
results_dir=fullfile(root,'results'); if ~exist(results_dir,'dir'); mkdir(results_dir); end
matfile=fullfile(results_dir,sprintf('exp7_algorithm_convergence_%s.mat',lower(RUN_MODE)));
save(matfile,'cfg','labels','main_idx','baseline_idx','proposed_idx','gamma_dB','gamma','P_total','P_total_dBm', ...
    'Obj','Pcomm','Feasible','Runtime','RuntimeMedian','RuntimeQ1','RuntimeQ3','HistRadar','HistMinPc','RadarProgress', ...
    'Stage1ExactEvals','Stage1Accepts','RefineExactEvals','TotalExactEvals','ProxyScreened', ...
    'SolverMode','RuntimeAccountingPolicy','SolverTelemetry','SolverInnerCalls','SolverDualitySuccess','SolverCVXAttempts', ...
    'SolverFallbacks','SolverFallbackRate','SolverMeanFPIIterations','SolverMeanPrimalDualGap', ...
    'SolverWallTime','InnerSolveMedian','CVXAttemptsMedian','FallbacksMedian','FallbackRateMedian', ...
    'FPIIterationsMedian','PrimalDualGapMedian','SolverWallMedian', ...
    'RadarGainTrace_dB','rad_gain_med','rad_gain_q1','rad_gain_q3','MinPcProgress','rad_updates','min_updates','CommonMain','ObjCommonMean_dB', ...
    'ObjCommonLo_dB','ObjCommonHi_dB','PcommCommonMean_dBm','PcommCommonLo_dBm', ...
    'PcommCommonHi_dBm','GainMean_dB','GainLo_dB','GainHi_dB','GainN', ...
    'Tsummary','Tpaired','Trefine','Truntime','RefineGain_dB','RefineAccepts', ...
    'Stage1Runtime','Stage2Runtime','BCSRefineRuntime','BCSExactEvals','V2Metrics','Tv2','mc_runtime');
csvfile=fullfile(results_dir,sprintf('exp7_algorithm_summary_%s.csv',lower(RUN_MODE)));
writetable(Tsummary,csvfile);
paired_csv=fullfile(results_dir,sprintf('exp7_paired_objective_gaps_%s.csv',lower(RUN_MODE)));
writetable(Tpaired,paired_csv);
refine_csv=fullfile(results_dir,sprintf('exp7_ma_rcg_mf_refinement_%s.csv',lower(RUN_MODE)));
writetable(Trefine,refine_csv);
runtime_csv=fullfile(results_dir,sprintf('exp7_runtime_breakdown_%s.csv',lower(RUN_MODE)));
writetable(Truntime,runtime_csv);
v2_csv=fullfile(results_dir,sprintf('exp7_v2_stage2_metrics_%s.csv',lower(RUN_MODE)));
writetable(Tv2,v2_csv);

fprintf('\nExperiment 7 completed in %.1f s (%.2f min).\n',mc_runtime,mc_runtime/60);
fprintf('Median accepted RadarSNR updates: %.1f\n',nanmedian_scalar(rad_updates));
fprintf('Median accepted MinPc updates: %.1f\n',nanmedian_scalar(min_updates));
fprintf('MAT: %s\nSummary CSV: %s\nPaired-gap CSV: %s\nRefinement CSV: %s\nRuntime CSV: %s\nV2 Stage-II CSV: %s\n', ...
    matfile,csvfile,paired_csv,refine_csv,runtime_csv,v2_csv);
if SAVE_FIGURES; fprintf('Figures: %s\n',figure_dir); else; fprintf('Per-experiment figures: not exported (%s profile)\n',FIG_POLICY.profile); end
fprintf('External RA code numerical baseline included: %d\n',cfg.exp7_compare_external_RA_code);
fprintf('Reason: %s\n',cfg.exp7_external_RA_reason);
fprintf('Total script runtime: %.1f s\n',toc(t_total));

%% ========================================================================
function [obj,pcomm,fea,trun,hRadar,hMinPc,refineGain,refineEvals,refineAccepts, ...
    stage1Evals,stage1Accepts,totalEvals,proxyScreened,stage1Runtime,stage2Runtime, ...
    bcsRefineRuntime,bcsExactEvals,v2metrics,solverMat]=run_one_exp7(cfg,mc,gamma,P_total,maxHist)
seed=cfg.seed_base+cfg.seed_stride*7+mc;
rng(seed,'twister');
scen=generate_isac_scenario(cfg,cfg.N,cfg.K,seed,[],[],struct());
obj=nan(1,8); pcomm=nan(1,8); fea=false(1,8); trun=nan(1,8);
refineGain=NaN; refineEvals=NaN; refineAccepts=NaN;
stage1Evals=NaN; stage1Accepts=NaN; totalEvals=NaN; proxyScreened=NaN;
stage1Runtime=NaN; stage2Runtime=NaN; bcsRefineRuntime=NaN; bcsExactEvals=NaN;
v2metrics=nan(1,8);
hRadar=nan(1,maxHist); hMinPc=nan(1,maxHist);
solverMat=nan(8,8);
comm_solver_telemetry('reset');

% Unified solutions. Use an explicit solver seed independent of channel draws.
rng(seed+700000,'twister');
t=tic; out=evaluate_ret_ris_schemes(cfg,scen,gamma,P_total); %#ok<NASGU>
obj(1)=out.Obj(1); pcomm(1)=out.Pcomm(1); fea(1)=out.Feasible(1); trun(1)=out.RuntimeCumulative(1);
obj(6)=out.Obj(4); pcomm(6)=out.Pcomm(4); fea(6)=out.Feasible(4); trun(6)=out.RuntimeCumulative(4);
if isfield(out,'PreRefineProposed') && out.PreRefineProposed.Feasible
    obj(7)=out.PreRefineProposed.Obj; pcomm(7)=out.PreRefineProposed.Pcomm;
    fea(7)=true; trun(7)=out.PreRefineProposed.RuntimeCumulative;
end
obj(8)=out.Obj(5); pcomm(8)=out.Pcomm(5); fea(8)=out.Feasible(5); trun(8)=out.RuntimeCumulative(5);
if isfield(out,'PreRefineProposed') && out.PreRefineProposed.Feasible
    stage1Runtime=out.PreRefineProposed.RuntimeCumulative;
end
if isfield(out,'ProposedMeta') && ~isempty(out.ProposedMeta)
    if isfield(out.ProposedMeta,'refinement_runtime_s'); stage2Runtime=out.ProposedMeta.refinement_runtime_s; end
    refineGain=out.ProposedMeta.refinement_gain_dB;
    refineEvals=out.ProposedMeta.total_exact_refine_evals;
    refineAccepts=double(out.ProposedMeta.mechanism_accepted)+out.ProposedMeta.polish_accepts;
    if isfield(out.ProposedMeta,'stage1_exact_evals'); stage1Evals=out.ProposedMeta.stage1_exact_evals; end
    if isfield(out.ProposedMeta,'stage1_accepts'); stage1Accepts=out.ProposedMeta.stage1_accepts; end
    if isfield(out.ProposedMeta,'total_exact_algorithm_evals'); totalEvals=out.ProposedMeta.total_exact_algorithm_evals; end
    if isfield(out.ProposedMeta,'stage1_screened') && isfield(out.ProposedMeta,'polish_screened')
        proxyScreened=out.ProposedMeta.stage1_screened+out.ProposedMeta.polish_screened;
    end
    v2metrics(1)=field_num(out.ProposedMeta,'v2_sweeps_completed');
    v2metrics(2)=field_num(out.ProposedMeta,'v2_block_visits');
    v2metrics(3)=field_num(out.ProposedMeta,'v2_budget_expansions');
    v2metrics(4)=field_num(out.ProposedMeta,'v2_ret_evals');
    v2metrics(5)=field_num(out.ProposedMeta,'v2_ris_evals');
    if isfield(out.ProposedMeta,'v2_trust_final_deg') && ~isempty(out.ProposedMeta.v2_trust_final_deg)
        z=out.ProposedMeta.v2_trust_final_deg; z=z(isfinite(z));
        if ~isempty(z); v2metrics(6)=mean(z); v2metrics(7)=max(z); end
    end
    if isfield(out.ProposedMeta,'v2_block_accept_count') && ~isempty(out.ProposedMeta.v2_block_accept_count)
        v2metrics(8)=sum(out.ProposedMeta.v2_block_accept_count);
    end
end

% Solver accounting for unified schemes uses the same cumulative path as runtime.
if isfield(out,'SolverStatsCumulative')
    solverMat(1,:)=comm_solver_stats_summary(out.SolverStatsCumulative{1});
    solverMat(6,:)=comm_solver_stats_summary(out.SolverStatsCumulative{4});
    if isfield(out,'PreRefineProposed') && isfield(out.PreRefineProposed,'SolverStatsCumulative')
        solverMat(7,:)=comm_solver_stats_summary(out.PreRefineProposed.SolverStatsCumulative);
    end
    solverMat(8,:)=comm_solver_stats_summary(out.SolverStatsCumulative{5});
end

% Sensing-oriented fixed-RET and sequential baselines with the full RIS candidate budget.
cfg_local=cfg;
cfg_local.exp6_fresh_candidate_mode='full';
cfg_local.exp6_fresh_allow_random=true;
cfg_local.exp6_fresh_ris_trials=cfg.num_phase_trials;
phi_seed=ones(cfg.N,1);

rng(seed+710001,'twister');
ss0=comm_solver_telemetry('snapshot');
t=tic; [~,r_fixed]=optimize_fixed_ret_ris_radar(cfg_local,scen,cfg.theta_ref, ...
    phi_seed,gamma,P_total); trun(2)=toc(t);
[obj(2),pcomm(2),fea(2)]=unpack_result(r_fixed);
ssx=comm_solver_telemetry('snapshot'); solverMat(2,:)=comm_solver_stats_summary(comm_solver_stats_delta(ss0,ssx));

theta_seq=out.Theta(3); if ~isfinite(theta_seq); theta_seq=cfg.theta_ref; end
rng(seed+710002,'twister');
ss0=comm_solver_telemetry('snapshot');
t=tic; [~,r_seq]=optimize_fixed_ret_ris_radar(cfg_local,scen,theta_seq, ...
    phi_seed,gamma,P_total); trun(3)=out.RuntimeCumulative(3)+toc(t);
[obj(3),pcomm(3),fea(3)]=unpack_result(r_seq);
ssx=comm_solver_telemetry('snapshot'); seqLocal=comm_solver_stats_delta(ss0,ssx);
if isfield(out,'SolverStatsBlock'); seqStats=comm_solver_stats_add(out.SolverStatsBlock{3},seqLocal); else; seqStats=seqLocal; end
solverMat(3,:)=comm_solver_stats_summary(seqStats);

% Random-search sanity baseline (supplementary only).
rng(seed+720000,'twister');
ss0=comm_solver_telemetry('snapshot');
t=tic; best=[]; best_obj=-Inf;
for q=1:cfg.exp7_random_trials
    ti=randi(numel(scen.theta_range)); theta_q=scen.theta_range(ti);
    phi_q=exp(1i*2*pi*rand(cfg.N,1));
    rr=evaluate_fixed_ret_ris_config(cfg,scen,theta_q,phi_q,gamma,P_total);
    if rr.feasible && rr.obj_snr>best_obj
        best=rr; best_obj=rr.obj_snr;
    end
end
trun(4)=toc(t);
ssx=comm_solver_telemetry('snapshot'); solverMat(4,:)=comm_solver_stats_summary(comm_solver_stats_delta(ss0,ssx));
if ~isempty(best); [obj(4),pcomm(4),fea(4)]=unpack_result(best); end

% Strong same-model derivative-free baseline.  To avoid weakening this
% comparator, it starts from the *same Stage-I radar-AO solution* used by the
% MA-RCG-MF refinement.  The runtime therefore includes that common
% Stage-I path plus the exhaustive block-coordinate refinement itself.
theta_bcs=cfg.theta_ref; phi_bcs=phi_seed; base_bcs_time=0;
if isfield(out,'PreRefineProposed') && out.PreRefineProposed.Feasible
    theta_bcs=out.PreRefineProposed.Theta; phi_bcs=out.PreRefineProposed.Phi;
    base_bcs_time=out.PreRefineProposed.RuntimeCumulative;
elseif out.Feasible(4)
    theta_bcs=out.Theta(4); phi_bcs=out.Phi{4}; base_bcs_time=out.RuntimeCumulative(4);
end
ss0=comm_solver_telemetry('snapshot');
t=tic; [r_bcs,~,bcsMeta]=optimize_ret_ris_block_coordinate(cfg,scen,gamma,P_total, ...
    theta_bcs,phi_bcs);
bcsRefineRuntime=toc(t); trun(5)=base_bcs_time+bcsRefineRuntime;
if isstruct(bcsMeta)&&isfield(bcsMeta,'exact_evals'); bcsExactEvals=bcsMeta.exact_evals; end
[obj(5),pcomm(5),fea(5)]=unpack_result(r_bcs);
ssx=comm_solver_telemetry('snapshot'); bcsLocal=comm_solver_stats_delta(ss0,ssx);
if isfield(out,'PreRefineProposed') && out.PreRefineProposed.Feasible && ...
        isfield(out.PreRefineProposed,'SolverStatsCumulative')
    bcsStats=comm_solver_stats_add(out.PreRefineProposed.SolverStatsCumulative,bcsLocal);
elseif isfield(out,'SolverStatsCumulative')
    bcsStats=comm_solver_stats_add(out.SolverStatsCumulative{4},bcsLocal);
else
    bcsStats=bcsLocal;
end
solverMat(5,:)=comm_solver_stats_summary(bcsStats);

if numel(out.History)>=5 && ~isempty(out.History{5})
    hRadar=pad_history(out.History{5}.ObjSNR,maxHist);
end
if numel(out.History)>=4 && ~isempty(out.History{4})
    hMinPc=pad_history(out.History{4}.Pcomm,maxHist);
end
end

function v=field_num(s,name)
v=NaN;
if isstruct(s) && isfield(s,name) && ~isempty(s.(name))
    z=s.(name);
    if isnumeric(z) && isscalar(z) && isfinite(z); v=double(z); end
end
end

function [o,p,f]=unpack_result(r)
f=isstruct(r)&&isfield(r,'feasible')&&r.feasible;
if f; o=r.obj_snr; p=r.p_comm; else; o=NaN; p=NaN; end
end

function y=pad_history(x,L)
y=nan(1,L); x=x(:).'; x=x(isfinite(x));
if isempty(x); return; end
n=min(numel(x),L); y(1:n)=x(1:n); y(n+1:end)=x(n);
end


function G=history_gain_db(H)
% Convert each true objective history to dB improvement over its own initial
% value. This shows accepted physical-objective gain without implying global
% optimality or normalizing every run to finish at one.
G=nan(size(H));
for ii=1:size(H,1)
    x=fill_terminal(H(ii,:));
    if isempty(x)||~isfinite(x(1))||x(1)<=0; continue; end
    ok=isfinite(x)&x>0;
    G(ii,ok)=10*log10(x(ok)/x(1));
end
end

function P=normalize_max_progress(H)
P=nan(size(H));
for ii=1:size(H,1)
    x=H(ii,:); ok=isfinite(x); if ~any(ok); continue; end
    x=fill_terminal(x); x0=x(1); xf=x(find(isfinite(x),1,'last'));
    den=xf-x0;
    if abs(den)<=1e-14; P(ii,isfinite(x))=1; else; P(ii,:)=(x-x0)/den; end
    P(ii,:)=min(max(P(ii,:),0),1);
end
end

function P=normalize_min_progress(H)
P=nan(size(H));
for ii=1:size(H,1)
    x=H(ii,:); ok=isfinite(x); if ~any(ok); continue; end
    x=fill_terminal(x); x0=x(1); xf=x(find(isfinite(x),1,'last'));
    den=x0-xf;
    if abs(den)<=1e-14; P(ii,isfinite(x))=1; else; P(ii,:)=(x0-x)/den; end
    P(ii,:)=min(max(P(ii,:),0),1);
end
end

function x=fill_terminal(x)
ids=find(isfinite(x)); if isempty(ids); return; end
x(ids(end)+1:end)=x(ids(end));
end

function n=count_effective_updates(H,is_max)
n=nan(size(H,1),1);
for ii=1:size(H,1)
    x=H(ii,:); x=x(isfinite(x)); if isempty(x); continue; end
    if is_max; d=diff(x)>1e-10*max(abs(x(1:end-1)),1e-12);
    else; d=diff(x)<-1e-10*max(abs(x(1:end-1)),1e-12); end
    n(ii)=1+sum(d);
end
end

function [med,q1,q3]=column_median_iqr(X)
med=nan(1,size(X,2)); q1=med; q3=med;
for jj=1:size(X,2)
    z=X(:,jj); z=z(isfinite(z)); if isempty(z); continue; end
    med(jj)=median(z); q1(jj)=empirical_quantile(z,0.25); q3(jj)=empirical_quantile(z,0.75);
end
end

function plot_band(x,med,q1,q3,marker)
valid=isfinite(med)&isfinite(q1)&isfinite(q3);
xv=x(valid); mv=med(valid); l=q1(valid); u=q3(valid);
fill([xv fliplr(xv)],[l fliplr(u)],[0.85 0.85 0.85], ...
    'EdgeColor','none','FaceAlpha',0.5); hold on;
plot(xv,mv,['-' marker],'LineWidth',1.8,'MarkerSize',6);
end

function [mu_db,lo_db,hi_db,n]=linear_mean_ci_to_db(X,mask,scale)
if nargin<3; scale=1; end
S=size(X,2); mu_db=nan(1,S); lo_db=nan(1,S); hi_db=nan(1,S); n=zeros(1,S);
for s=1:S
    x=X(:,s); ok=mask(:,s)&isfinite(x)&x>=0; x=x(ok)*scale; n(s)=numel(x);
    if isempty(x); continue; end
    mu=mean(x); se=std(x,0)/sqrt(max(n(s),1));
    if n(s)>1; d=student_t_critical(0.95,n(s)-1)*se; else; d=0; end
    lo=max(mu-d,realmin); hi=max(mu+d,realmin);
    mu_db(s)=10*log10(max(mu,realmin)); lo_db(s)=10*log10(lo); hi_db(s)=10*log10(hi);
end
end

function [mu,lo,hi,n]=paired_db_gain(Obj,Feasible,proposed_idx,baseline_idx,common_mask)
B=numel(baseline_idx); mu=nan(1,B); lo=mu; hi=mu; n=zeros(1,B);
if nargin<5 || isempty(common_mask); common_mask=true(size(Obj,1),1); end
for bb=1:B
    b=baseline_idx(bb);
    ok=common_mask & Feasible(:,proposed_idx)&Feasible(:,b)& ...
       isfinite(Obj(:,proposed_idx))&isfinite(Obj(:,b))& ...
       Obj(:,proposed_idx)>0&Obj(:,b)>0;
    d=10*log10(Obj(ok,proposed_idx)./Obj(ok,b));
    n(bb)=numel(d);
    if isempty(d); continue; end
    mu(bb)=mean(d);
    if numel(d)>1; half=student_t_critical(0.95,numel(d)-1)*std(d,0)/sqrt(numel(d)); else; half=0; end
    lo(bb)=mu(bb)-half; hi(bb)=mu(bb)+half;
end
end

function m=nanmedian_scalar(x)
x=x(isfinite(x));
if isempty(x); m=NaN; else; m=median(x); end
end

function m=nanmedian_local(X,dim)
if nargin<2; dim=1; end
if dim~=1; error('nanmedian_local currently supports dim=1 only.'); end
m=nan(1,size(X,2));
for j=1:size(X,2); z=X(:,j); z=z(isfinite(z)); if ~isempty(z); m(j)=median(z); end; end
end

function q=nanquantile_local(X,p,dim)
if nargin<3; dim=1; end
if dim~=1; error('nanquantile_local currently supports dim=1 only.'); end
q=nan(1,size(X,2));
for j=1:size(X,2); z=X(:,j); z=z(isfinite(z)); if ~isempty(z); q(j)=empirical_quantile(z,p); end; end
end

function save_exp7_figure(fig,folder,name)
pngfile=fullfile(folder,[name '.png']); figfile=fullfile(folder,[name '.fig']);
try; exportgraphics(fig,pngfile,'Resolution',300); catch; print(fig,pngfile,'-dpng','-r300'); end
savefig(fig,figfile);
end

function exp7_progress(total)
persistent count t0
if isempty(count)||count>=total; count=0; t0=tic; end
count=count+1; elapsed=toc(t0); eta=elapsed/max(count,1)*(total-count);
fprintf('\rExp7 progress: %3.0f%% (%d/%d), elapsed %.1f min, ETA %.1f min', ...
    100*count/total,count,total,elapsed/60,eta/60);
if count>=total; fprintf('\n'); end
end
