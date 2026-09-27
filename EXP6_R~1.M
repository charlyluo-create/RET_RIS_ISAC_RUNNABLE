% =========================================================================
% Experiment 6 FINAL: Multi-Timescale Mobility and Configuration Staleness
%
% Purpose
%   Quantify how vehicle speed and RIS update period affect the performance
%   of a roadside RET-RIS ISAC system under a practical multi-timescale
%   control architecture.
%
% Control timescales
%   Slow   : RET tilt is selected for the traffic/service state and held.
%   Medium : RIS phase is refreshed periodically every tau_RIS milliseconds.
%   Fast   : communication and sensing beamformers are re-optimized from
%            the current CSI for every evaluated state.
%
% Comparison at t0 + tau_RIS
%   Stale RIS : keep the RET/RIS configuration selected at t0; update Wc/Wr.
%   Fresh RIS : keep the same slow RET tilt, refresh RIS using the current
%               geometry/channel, and update Wc/Wr.
%
% Main metrics
%   Mean and 90th-percentile Delta_SNR(v,tau), where
%   Delta_SNR = SNR_fresh_RIS - SNR_stale_RIS [dB].
%
% Strict paired controls
%   One initial geometry, motion-direction vector, target-motion vector and
%   normalized fading realization are drawn per MC replicate and reused for
%   every speed/update-period cell (common random numbers).  All moved channel
%   states are precomputed before any optimizer call.  tau_RIS = 0 ms is an
%   exact calibration point.  For tau_RIS > 0 the fresh RIS benchmark is
%   deterministic, warm-started from the stale phase and contains no random
%   restart.  P90 confidence intervals use one shared bootstrap index vector
%   per replicate so adjacent-period comparisons remain paired.
%
% This experiment isolates geometry-induced passive-configuration staleness.
% It does not claim symbol-level tracking of Doppler-varying small-scale CSI;
% the active beamformers are assumed to handle the fast CSI timescale.
% =========================================================================

clear; close all; clc;

%% ======================== Controls =======================================
env_mode=getenv('RET_RIS_RUN_MODE');
if isempty(env_mode); RUN_MODE='paper'; else; RUN_MODE=lower(env_mode); end
SHOW_PROGRESS=true;
t_total=tic;

cfg=unified_ret_ris_isac_config(RUN_MODE);
FIG_POLICY=get_figure_policy(RUN_MODE);
SAVE_MAIN_FIGURE=FIG_POLICY.save_main;
SAVE_SUPP_FIGURES=FIG_POLICY.save_supplementary;
SAVE_FIGURES=FIG_POLICY.save_any;
ITER=cfg.ITER_EXP6;
speed_kmh=cfg.speed_kmh_exp6;
update_ms=cfg.ris_update_ms_exp6;
numV=numel(speed_kmh);
numTau=numel(update_ms);
gamma_dB=cfg.gamma_exp6_dB;
gamma=10^(gamma_dB/10);
P_total_dBm=cfg.P_total_exp6_dBm;
P_total=cfg.P_total_exp6;
num_workers=min(cfg.num_workers,ITER);

%% ======================== Paths / Checks =================================
root=fileparts(mfilename('fullpath'));
if isempty(root); root=pwd; end
addpath(root);
func_path=fullfile(root,'function');
if exist(func_path,'dir'); addpath(func_path); end

required={'generate_isac_scenario','evaluate_ret_ris_schemes', ...
    'evaluate_fixed_ret_ris_config','optimize_fixed_ret_ris_radar', ...
    'generate_channel','solve_comm_min_power','build_radar_matrix', ...
    'build_effective_user_channel','allocate_radar_power','apply_relative_motion_phase', ...
    'bootstrap_quantile_grid_paired','paired_adjacent_trend_stats'};
for ii=1:numel(required)
    if exist(required{ii},'file')~=2
        error('Required function not found: %s.m',required{ii});
    end
end
if exist('cvx_begin','file')~=2
    error('CVX is required. Install CVX and run cvx_setup first.');
end

stamp=datestr(now,'yyyymmdd_HHMMSS');
figure_dir=fullfile(root,sprintf('Exp6_Figures_%s_%s',lower(RUN_MODE),stamp));
if SAVE_FIGURES && ~exist(figure_dir,'dir'); mkdir(figure_dir); end

%% ======================== Parallel Pool ==================================
use_parallel=true;
if use_parallel
    try
        poolobj=gcp('nocreate');
        if ~isempty(poolobj)&&poolobj.NumWorkers~=num_workers
            delete(poolobj); poolobj=[];
        end
        if isempty(poolobj); poolobj=parpool(num_workers); end
        pctRunOnAll(['addpath(''' root ''');']);
        pctRunOnAll(['addpath(''' func_path ''');']);
    catch ME
        warning('Parallel setup failed: %s. Falling back to serial.',ME.message);
        use_parallel=false;
    end
end

max_disp=max(speed_kmh)/3.6*max(update_ms)/1000;
fprintf('\n=================================================================\n');
fprintf(' EXPERIMENT 6: MULTI-TIMESCALE MOBILITY / RIS STALENESS\n');
fprintf(' Mode=%s, MC=%d, gamma=%.1f dB, P_T=%.1f dBm\n', ...
    upper(RUN_MODE),ITER,gamma_dB,P_total_dBm);
fprintf(' Vehicle speeds (km/h): '); fprintf('%g ',speed_kmh); fprintf('\n');
fprintf(' RIS update periods (ms): '); fprintf('%g ',update_ms); fprintf('\n');
fprintf(' Maximum displacement in one update interval: %.3f m\n',max_disp);
fprintf(' Slow RET held fixed; fresh benchmark refreshes RIS only.\n');
fprintf(' Communication-outage stress point: gamma=%.1f dB, P_T=%.1f dBm\n', ...
    cfg.exp6_stress_gamma_dB,cfg.exp6_stress_P_total_dBm);
if SAVE_FIGURES; fprintf(' Figures: %s\n',figure_dir); end
fprintf('=================================================================\n\n');

%% ======================== Buffers ========================================
Loss_dB=nan(ITER,numV,numTau);
RawLoss_dB=nan(ITER,numV,numTau);
ObjStale=nan(ITER,numV,numTau);
ObjFresh=nan(ITER,numV,numTau);
PcommStale=nan(ITER,numV,numTau);
PcommFresh=nan(ITER,numV,numTau);
CommChange_dB=nan(ITER,numV,numTau);
PhiDistance=nan(ITER,numV,numTau);
StaleFeasible=false(ITER,numV,numTau);
FreshFeasible=false(ITER,numV,numTau);
StressStaleFeasible=false(ITER,numV,numTau);
StressFreshFeasible=false(ITER,numV,numTau);
FreshDominanceFallback=false(ITER,numV,numTau);
InitialFeasible=false(ITER,1);
InitialTheta=nan(ITER,1);
% Pairing metadata: one traffic/channel realization per MC is shared by all
% speed/update-period cells.  These arrays make the pairing auditable.
InitialX=nan(ITER,cfg.K); InitialY=nan(ITER,cfg.K);
UserMotionDirection=nan(ITER,cfg.K); TargetMotionDirection=nan(ITER,cfg.T);
ChannelSeed=nan(ITER,1); PairingRebuildError=nan(ITER,1);
CompletedMC=false(ITER,1);

%% ======================== Progress =======================================
progressQueue=[];
if SHOW_PROGRESS
    try
        progressQueue=parallel.pool.DataQueue;
        afterEach(progressQueue,@exp6_progress_bar);
    catch ME
        warning('Progress bar initialization failed: %s',ME.message);
        progressQueue=[];
    end
end

%% ======================== Monte Carlo ====================================
% Exp.6 is deliberately run in checkpointed batches. PAPER mode uses more
% MC samples than earlier releases because tail statistics (P90) are noisy.
% Batch checkpointing makes the larger run restart-safe without changing any
% sample identity or random seed.
results_dir=fullfile(root,'results');
if ~exist(results_dir,'dir'); mkdir(results_dir); end
checkpoint_file=fullfile(results_dir,sprintf('exp6_multitimescale_mobility_%s_checkpoint.mat',lower(RUN_MODE)));
CheckpointSignature=exp6_checkpoint_signature(cfg,speed_kmh,update_ms,gamma,P_total); %#ok<NASGU>
if cfg.exp6_resume_checkpoint && exist(checkpoint_file,'file')
    try
        C=load(checkpoint_file,'CP');
        if isfield(C,'CP') && isfield(C.CP,'signature') && strcmp(C.CP.signature,CheckpointSignature)
            CP=C.CP;
            names={'Loss_dB','RawLoss_dB','ObjStale','ObjFresh','PcommStale','PcommFresh', ...
                'CommChange_dB','PhiDistance','StaleFeasible','FreshFeasible', ...
                'StressStaleFeasible','StressFreshFeasible','FreshDominanceFallback', ...
                'InitialFeasible','InitialTheta','InitialX','InitialY','UserMotionDirection', ...
                'TargetMotionDirection','ChannelSeed','PairingRebuildError','CompletedMC'};
            for jj=1:numel(names)
                if isfield(CP,names{jj}); eval([names{jj} '=CP.(names{jj});']); end %#ok<EVLDIR>
            end
            fprintf('Resuming Exp.6 checkpoint: %d/%d MC replicates already complete.\n',nnz(CompletedMC),ITER);
        else
            warning('Ignoring Exp.6 checkpoint because its configuration signature does not match.');
        end
    catch ME
        warning('Could not resume Exp.6 checkpoint: %s. Starting from the available in-memory buffers.',ME.message);
    end
end

pending=find(~CompletedMC);
if SHOW_PROGRESS && ~isempty(progressQueue)
    exp6_progress_bar(struct('action','reset','total',max(1,numel(pending)*numV*numTau)));
end
t_mc=tic;
batch_size=max(1,min(cfg.exp6_checkpoint_batch,max(1,numel(pending))));
for b0=1:batch_size:numel(pending)
    ids=pending(b0:min(b0+batch_size-1,numel(pending)));
    batchOut=cell(numel(ids),1);
    if use_parallel
        parfor jj=1:numel(ids)
            iter=ids(jj); seed=cfg.seed_base+iter*cfg.seed_stride;
            batchOut{jj}=run_one_mc_mobility(cfg,speed_kmh,update_ms,gamma,P_total, ...
                seed,max_disp,progressQueue);
        end
    else
        for jj=1:numel(ids)
            iter=ids(jj); seed=cfg.seed_base+iter*cfg.seed_stride;
            batchOut{jj}=run_one_mc_mobility(cfg,speed_kmh,update_ms,gamma,P_total, ...
                seed,max_disp,progressQueue);
        end
    end
    for jj=1:numel(ids)
        iter=ids(jj); o=batchOut{jj};
        Loss_dB(iter,:,:)=reshape(o.Loss_dB,[1 numV numTau]);
        RawLoss_dB(iter,:,:)=reshape(o.RawLoss_dB,[1 numV numTau]);
        ObjStale(iter,:,:)=reshape(o.ObjStale,[1 numV numTau]);
        ObjFresh(iter,:,:)=reshape(o.ObjFresh,[1 numV numTau]);
        PcommStale(iter,:,:)=reshape(o.PcommStale,[1 numV numTau]);
        PcommFresh(iter,:,:)=reshape(o.PcommFresh,[1 numV numTau]);
        CommChange_dB(iter,:,:)=reshape(o.CommChange_dB,[1 numV numTau]);
        PhiDistance(iter,:,:)=reshape(o.PhiDistance,[1 numV numTau]);
        StaleFeasible(iter,:,:)=reshape(o.StaleFeasible,[1 numV numTau]);
        FreshFeasible(iter,:,:)=reshape(o.FreshFeasible,[1 numV numTau]);
        StressStaleFeasible(iter,:,:)=reshape(o.StressStaleFeasible,[1 numV numTau]);
        StressFreshFeasible(iter,:,:)=reshape(o.StressFreshFeasible,[1 numV numTau]);
        FreshDominanceFallback(iter,:,:)=reshape(o.FreshDominanceFallback,[1 numV numTau]);
        InitialFeasible(iter)=o.InitialFeasible; InitialTheta(iter)=o.InitialTheta;
        InitialX(iter,:)=o.InitialX; InitialY(iter,:)=o.InitialY;
        UserMotionDirection(iter,:)=o.UserMotionDirection;
        TargetMotionDirection(iter,:)=o.TargetMotionDirection;
        ChannelSeed(iter)=o.ChannelSeed; PairingRebuildError(iter)=o.PairingRebuildError;
        CompletedMC(iter)=true;
    end
    if cfg.exp6_resume_checkpoint
        CP=exp6_checkpoint_struct(CheckpointSignature,Loss_dB,RawLoss_dB,ObjStale,ObjFresh, ...
            PcommStale,PcommFresh,CommChange_dB,PhiDistance,StaleFeasible,FreshFeasible, ...
            StressStaleFeasible,StressFreshFeasible,FreshDominanceFallback,InitialFeasible, ...
            InitialTheta,InitialX,InitialY,UserMotionDirection,TargetMotionDirection, ...
            ChannelSeed,PairingRebuildError,CompletedMC); %#ok<NASGU>
        save(checkpoint_file,'CP','-v7.3');
        fprintf('Exp.6 checkpoint saved: %d/%d MC complete.\n',nnz(CompletedMC),ITER);
    end
end
runtime=toc(t_mc);
%% ======================== Statistics =====================================
[LossMean,LossCI,LossCount]=mean_ci_dim1(Loss_dB);
if cfg.exp6_use_paired_bootstrap
    [LossP90,LossP90Lo,LossP90Hi,LossP90Count,LossP90Boot]=bootstrap_quantile_grid_paired( ...
        Loss_dB,0.90,cfg.exp6_bootstrap_reps,cfg.exp6_bootstrap_seed);
else
    [LossP90,LossP90Lo,LossP90Hi,LossP90Count]=bootstrap_quantile_grid( ...
        Loss_dB,0.90,cfg.exp6_bootstrap_reps,cfg.exp6_bootstrap_seed);
    LossP90Boot=[];
end
if cfg.exp6_trend_audit
    TrendAudit=paired_adjacent_trend_stats(Loss_dB,LossP90,LossP90Boot,update_ms);
else
    TrendAudit=struct();
end
[CommChangeMean,CommChangeCI,CommCount]=mean_ci_dim1(CommChange_dB);
[PhiDistanceMean,~,PhiDistanceCount]=mean_ci_dim1(PhiDistance);
zero_idx=find(update_ms==0,1);
if isempty(zero_idx)
    ZeroCalibrationMaxAbs_dB=NaN;
else
    z=RawLoss_dB(:,:,zero_idx); z=z(isfinite(z));
    if isempty(z); ZeroCalibrationMaxAbs_dB=NaN; else; ZeroCalibrationMaxAbs_dB=max(abs(z)); end
end
StaleRate=squeeze(100*mean(StaleFeasible,1));
FreshRate=squeeze(100*mean(FreshFeasible,1));

% Conditional outage is reported at both the nominal operating point and an
% explicitly labelled communication-stress point. It excludes states that are
% infeasible even after a fresh RIS update.
StaleOutage=conditional_outage(FreshFeasible,StaleFeasible);
StressOutage=conditional_outage(StressFreshFeasible,StressStaleFeasible);

if SAVE_MAIN_FIGURE
%% ======================== Main Figure ====================================
fig1=figure('Position',[40 80 1450 570],'Color','w');
tl=tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile;
imagesc(1:numTau,1:numV,LossMean);
axis xy; box on;
xticks(1:numTau); xticklabels(string(update_ms));
yticks(1:numV); yticklabels(string(speed_kmh));
xlabel('RIS update period \tau_{RIS} (ms)');
ylabel('Vehicle speed v (km/h)');
title('(a) Mean paired radar-SNR staleness loss');
cb=colorbar; cb.Label.String='Mean radar-SNR loss (dB)';
add_heatmap_values(1:numTau,1:numV,LossMean,'%.2f');

nexttile;
imagesc(1:numTau,1:numV,LossP90);
axis xy; box on;
xticks(1:numTau); xticklabels(string(update_ms));
yticks(1:numV); yticklabels(string(speed_kmh));
xlabel('RIS update period \tau_{RIS} (ms)');
ylabel('Vehicle speed v (km/h)');
title(sprintf('(b) 90th-percentile paired loss (%d paired-bootstrap replicates)', ...
    cfg.exp6_bootstrap_reps));
cb=colorbar; cb.Label.String='90th-percentile radar-SNR loss (dB)';
caxis([0 max(0.05,finite_max(LossP90))]);
add_heatmap_values(1:numTau,1:numV,LossP90,'%.2f');

tt=title(tl,sprintf(['Strict-Paired Multi-Timescale RET-RIS Aging: Vehicle Mobility vs. ', ...
    'RIS Update Period (gamma = %.1f dB, P_T = %.1f dBm)'], ...
    gamma_dB,P_total_dBm));
tt.Interpreter='none';
if SAVE_MAIN_FIGURE
    save_figure_auto(fig1,figure_dir,'01_Vehicle_Speed_Update_Period_Performance_Loss');
end

%% ======================== Supplementary Curves ===========================
end

if SAVE_SUPP_FIGURES
fig2=figure('Position',[80 80 980 630],'Color','w');
markers={'-o','-s','-d','-^','-v','-p'};
for vi=1:numV
    mk=markers{1+mod(vi-1,numel(markers))};
    errorbar(update_ms,LossMean(vi,:),LossCI(vi,:),mk, ...
        'LineWidth',1.6,'MarkerSize',8); hold on;
end
xlabel('RIS update period \tau_{RIS} (ms)');
ylabel('Mean radar-SNR loss (dB)');
legend(compose('v = %g km/h',speed_kmh),'Location','northwest');
grid on;
title('Passive-Configuration Aging under Vehicle Mobility');
if SAVE_SUPP_FIGURES
    save_figure_auto(fig2,figure_dir,'S1_Radar_SNR_Loss_vs_Update_Period');
end

fig3=figure('Position',[100 100 980 630],'Color','w');
imagesc(1:numTau,1:numV,CommChangeMean); axis xy; box on;
xticks(1:numTau); xticklabels(string(update_ms));
yticks(1:numV); yticklabels(string(speed_kmh));
xlabel('RIS update period \tau_{RIS} (ms)'); ylabel('Vehicle speed v (km/h)');
title('Supplementary Communication-Power Change with Stale RIS');
cb=colorbar; cb.Label.String='Mean 10log_{10}(P_c^{stale}/P_c^{fresh}) (dB)';
clim=max(0.05,finite_abs_max(CommChangeMean)); caxis([-clim clim]);
add_heatmap_values(1:numTau,1:numV,CommChangeMean,'%.2f');
if SAVE_SUPP_FIGURES
    save_figure_auto(fig3,figure_dir,'S2_Communication_Power_Change');
end

fig4=figure('Position',[120 120 980 630],'Color','w');
imagesc(1:numTau,1:numV,StressOutage); axis xy; box on;
xticks(1:numTau); xticklabels(string(update_ms));
yticks(1:numV); yticklabels(string(speed_kmh));
xlabel('RIS update period \tau_{RIS} (ms)'); ylabel('Vehicle speed v (km/h)');
title(sprintf('Supplementary Communication Outage at Stress Point (QoS = %.0f dB, P_T = %.0f dBm)', ...
    cfg.exp6_stress_gamma_dB,cfg.exp6_stress_P_total_dBm));
cb=colorbar; cb.Label.String='Conditional outage rate (%)';
caxis([0 max(1,finite_max(StressOutage))]);
add_heatmap_values(1:numTau,1:numV,StressOutage,'%.1f');
if SAVE_SUPP_FIGURES
    save_figure_auto(fig4,figure_dir,'S3_Stress_Point_Communication_Outage');
end


%% ======================== Displacement-Domain Diagnostic =================
if cfg.exp6_report_displacement_curve
    fig5=figure('Position',[140 140 980 630],'Color','w'); hold on;
    markers={'-o','-s','-d','-^','-v','-p'};
    for vi=1:numV
        displacement_m=(speed_kmh(vi)/3.6)*(update_ms/1000);
        errorbar(displacement_m,LossMean(vi,:),LossCI(vi,:), ...
            markers{1+mod(vi-1,numel(markers))},'LineWidth',1.6,'MarkerSize',7);
    end
    grid on; xlabel('Vehicle displacement during one RIS interval v\tau_{RIS} (m)');
    ylabel('Mean radar-SNR loss (dB)');
    legend(compose('v = %g km/h',speed_kmh),'Location','northwest');
    title('Geometry-Induced RIS Aging versus Displacement');
    if SAVE_SUPP_FIGURES
        save_figure_auto(fig5,figure_dir,'S4_Radar_Loss_vs_Displacement');
    end

    fig6=figure('Position',[160 160 980 630],'Color','w');
    imagesc(1:numTau,1:numV,LossP90Hi-LossP90Lo); axis xy; box on;
    xticks(1:numTau); xticklabels(string(update_ms));
    yticks(1:numV); yticklabels(string(speed_kmh));
    xlabel('RIS update period \tau_{RIS} (ms)'); ylabel('Vehicle speed v (km/h)');
    title('Bootstrap 95% CI Width of the 90th-Percentile Loss');
    cb=colorbar; cb.Label.String='CI width (dB)';
    add_heatmap_values(1:numTau,1:numV,LossP90Hi-LossP90Lo,'%.2f');
    if SAVE_SUPP_FIGURES
        save_figure_auto(fig6,figure_dir,'S5_P90_Bootstrap_CI_Width');
    end
end

if cfg.exp6_trend_audit && isfield(TrendAudit,'MeanDelta')
    fig7=figure('Position',[180 180 1120 700],'Color','w');
    tl7=tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
    xt=1:(numTau-1); pairLabels=compose('%g->%g',update_ms(1:end-1),update_ms(2:end));
    nexttile; hold on;
    for vi=1:numV
        errorbar(xt,TrendAudit.MeanDelta(vi,:), ...
            TrendAudit.MeanDelta(vi,:)-TrendAudit.MeanDeltaLo(vi,:), ...
            TrendAudit.MeanDeltaHi(vi,:)-TrendAudit.MeanDelta(vi,:),'-o','LineWidth',1.2);
    end
    yline(0,'--k'); grid on; set(gca,'XTick',xt,'XTickLabel',pairLabels);
    ylabel('Paired mean-loss change (dB)'); title('(a) Adjacent update-period mean changes');
    legend(compose('v=%g km/h',speed_kmh),'Location','best');
    nexttile; hold on;
    for vi=1:numV
        errorbar(xt,TrendAudit.QuantileDelta(vi,:), ...
            TrendAudit.QuantileDelta(vi,:)-TrendAudit.QuantileDeltaLo(vi,:), ...
            TrendAudit.QuantileDeltaHi(vi,:)-TrendAudit.QuantileDelta(vi,:),'-s','LineWidth',1.2);
    end
    yline(0,'--k'); grid on; set(gca,'XTick',xt,'XTickLabel',pairLabels);
    xlabel('Adjacent sampled update periods (ms)'); ylabel('Paired P90 change (dB)');
    title('(b) Adjacent update-period P90 changes (paired bootstrap)');
    title(tl7,'Paired Trend Audit; Negative Values Are Not Removed or Smoothed');
    if SAVE_SUPP_FIGURES
        save_figure_auto(fig7,figure_dir,'S6_Paired_Adjacent_Trend_Audit');
    end
end

end

%% ======================== Console / Save =================================
fprintf('\nMean radar-SNR loss (dB), rows=speed, columns=update period:\n');
disp(array2table([speed_kmh(:),LossMean], ...
    'VariableNames',[{'Speed_kmh'},make_valid_names(update_ms,'Tau_ms_')]));
fprintf('\nMean communication-power change 10log10(Pc_stale/Pc_fresh) (dB):\n');
disp(array2table([speed_kmh(:),CommChangeMean], ...
    'VariableNames',[{'Speed_kmh'},make_valid_names(update_ms,'Tau_ms_')]));

fprintf('\nNominal staleness-induced outage (%%):\n');
disp(array2table([speed_kmh(:),StaleOutage], ...
    'VariableNames',[{'Speed_kmh'},make_valid_names(update_ms,'Tau_ms_')]));
fprintf('\nStress-point staleness-induced outage (%%):\n');
disp(array2table([speed_kmh(:),StressOutage], ...
    'VariableNames',[{'Speed_kmh'},make_valid_names(update_ms,'Tau_ms_')]));
fprintf('Initial RET-RIS-RadarSNR feasibility: %.1f %%\n',100*mean(InitialFeasible));
fprintf('Maximum mean global-phase-invariant RIS refresh distance: %.4f\n', ...
    finite_max(PhiDistanceMean));
fprintf('Fresh-dominance safeguard activations: %d\n',nnz(FreshDominanceFallback));
fprintf('Minimum paired feasible sample count in the loss grid: %d\n',min(LossCount(:)));

fprintf('Minimum bootstrap sample count for P90: %d\n',min(LossP90Count(:)));
fprintf('Maximum bootstrap 95%% CI width for P90: %.3f dB\n', ...
    finite_max(LossP90Hi-LossP90Lo));
fprintf('Pairing mode: %s\n',cfg.exp6_pairing_version);
fprintf('Maximum state-bank rebuild error: %.3e\n',finite_max(PairingRebuildError));
if cfg.exp6_trend_audit && isfield(TrendAudit,'NumPointMeanDecreases')
    fprintf(['Adjacent-period mean decreases: %d point estimates; %d statistically ', ...
        'resolved at paired 95%% CI.\n'],TrendAudit.NumPointMeanDecreases, ...
        TrendAudit.NumSignificantMeanDecreases);
    fprintf(['Adjacent-period P90 decreases: %d point estimates; %d statistically ', ...
        'resolved by paired bootstrap 95%% CI.\n'],TrendAudit.NumPointQuantileDecreases, ...
        TrendAudit.NumSignificantQuantileDecreases);
end

if any(update_ms>0) && finite_max(LossMean(:,update_ms>0))<1e-6
    warning(['All positive-period sensing losses are numerically zero. ', ...
        'Inspect target motion, path-phase evolution and fresh-RIS acceptance diagnostics.']);
end
if ~isempty(zero_idx)
    fprintf('Zero-ms calibration max |raw loss|: %.3e dB\n',ZeroCalibrationMaxAbs_dB);
end

results_dir=fullfile(root,'results');
if ~exist(results_dir,'dir'); mkdir(results_dir); end
% Export an explicit paired trend table. Negative adjacent changes are retained;
% the CI columns tell whether a local reversal is statistically resolved.
Ttrend=trend_audit_table(TrendAudit,speed_kmh,update_ms);
writetable(Ttrend,fullfile(results_dir,sprintf('exp6_paired_trend_audit_%s.csv',lower(RUN_MODE))));
audit_file=fullfile(results_dir,sprintf('exp6_pairing_audit_%s.txt',lower(RUN_MODE)));
write_pairing_audit(audit_file,cfg,ITER,speed_kmh,update_ms,LossCount,LossP90Count, ...
    PairingRebuildError,FreshDominanceFallback,TrendAudit,ZeroCalibrationMaxAbs_dB);

save(fullfile(results_dir,sprintf('exp6_multitimescale_mobility_%s.mat',lower(RUN_MODE))), ...
    'cfg','speed_kmh','update_ms','gamma_dB','gamma','P_total','P_total_dBm', ...
    'Loss_dB','RawLoss_dB','CommChange_dB','PhiDistance','ObjStale','ObjFresh','PcommStale','PcommFresh', ...
    'StaleFeasible','FreshFeasible','StressStaleFeasible','StressFreshFeasible','FreshDominanceFallback', ...
    'InitialFeasible','InitialTheta','InitialX','InitialY','UserMotionDirection','TargetMotionDirection', ...
    'ChannelSeed','PairingRebuildError','CompletedMC','LossMean','LossCI','LossP90','LossP90Lo', ...
    'LossP90Hi','LossP90Count','LossCount','TrendAudit','Ttrend', ...
    'CommChangeMean','CommChangeCI','CommCount','PhiDistanceMean','PhiDistanceCount', ...
    'ZeroCalibrationMaxAbs_dB','StaleRate','FreshRate','StaleOutage','StressOutage','runtime');
if cfg.exp6_delete_checkpoint_on_success && exist(checkpoint_file,'file')
    delete(checkpoint_file);
end

total_time=toc(t_total);
fprintf('\n=================================================================\n');
fprintf(' EXPERIMENT 6 COMPLETED\n');
fprintf(' Monte Carlo time: %.1f s (%.2f min / %.2f h)\n', ...
    runtime,runtime/60,runtime/3600);
fprintf(' Total runtime:     %.1f s (%.2f min / %.2f h)\n', ...
    total_time,total_time/60,total_time/3600);
if SAVE_FIGURES; fprintf(' Figure folder:     %s\n',figure_dir); end
fprintf('=================================================================\n');


%% ========================================================================
%% Local functions
%% ========================================================================
function o=run_one_mc_mobility(cfg,speed_kmh,update_ms,gamma,P_total, ...
    seed,max_disp,progressQueue)

numV=numel(speed_kmh); numTau=numel(update_ms);
o.Loss_dB=nan(numV,numTau); o.RawLoss_dB=nan(numV,numTau);
o.ObjStale=nan(numV,numTau); o.ObjFresh=nan(numV,numTau);
o.PcommStale=nan(numV,numTau); o.PcommFresh=nan(numV,numTau);
o.CommChange_dB=nan(numV,numTau); o.PhiDistance=nan(numV,numTau);
o.StaleFeasible=false(numV,numTau); o.FreshFeasible=false(numV,numTau);
o.StressStaleFeasible=false(numV,numTau); o.StressFreshFeasible=false(numV,numTau);
o.FreshDominanceFallback=false(numV,numTau);
o.InitialFeasible=false; o.InitialTheta=NaN;
o.InitialX=nan(1,cfg.K); o.InitialY=nan(1,cfg.K);
o.UserMotionDirection=nan(1,cfg.K); o.TargetMotionDirection=nan(1,cfg.T);
o.ChannelSeed=NaN; o.PairingRebuildError=NaN;

if isfield(cfg,'exp6_antithetic_motion') && cfg.exp6_antithetic_motion
    error(['exp6_antithetic_motion is intentionally disabled in the locked manuscript profile. ', ...
        'The reported MC unit is one physical traffic realization.']);
end
if isfield(cfg,'exp6_enforce_deterministic_fresh') && cfg.exp6_enforce_deterministic_fresh && ...
        isfield(cfg,'exp6_fresh_allow_random') && cfg.exp6_fresh_allow_random
    error('Strict paired Exp.6 forbids random fresh-RIS restarts. Set exp6_fresh_allow_random=false.');
end

% One physical traffic/channel realization is drawn per MC.  Every cell in
% the speed/update-period grid reuses these latent variables.  The moved
% state bank is built before any optimization call so optimizer internals can
% never contaminate the channel random stream of another cell.
rng(seed,'twister');
x_lo=cfg.road_x_min+max_disp; x_hi=cfg.road_x_max-max_disp;
if x_hi<=x_lo; error('Road interval is too short for the selected speed/update grid.'); end
x0=x_lo+(x_hi-x_lo)*rand(1,cfg.K);
y0=-cfg.road_y_half+2*cfg.road_y_half*rand(1,cfg.K);
direction=2*(rand(1,cfg.K)>0.5)-1; direction(direction==0)=1;
target_direction=2*(rand(1,cfg.T)>0.5)-1; target_direction(target_direction==0)=1;
channel_seed=seed+100;
o.InitialX=x0; o.InitialY=y0; o.UserMotionDirection=direction;
o.TargetMotionDirection=target_direction; o.ChannelSeed=channel_seed;

scen0=generate_isac_scenario(cfg,cfg.N,cfg.K,channel_seed,x0,y0,struct());
stateBank=cell(numV,numTau);
for vi=1:numV
    v_mps=speed_kmh(vi)/3.6;
    for ti=1:numTau
        tau_s=update_ms(ti)/1000;
        if update_ms(ti)==0
            stateBank{vi,ti}=scen0;
        else
            stateBank{vi,ti}=build_paired_moved_snapshot(cfg,scen0,channel_seed,x0,y0, ...
                direction,target_direction,v_mps,tau_s);
        end
    end
end

% Rebuild one endpoint independently and verify bit-level/common-random-number
% reproducibility up to floating-point roundoff.
if isfield(cfg,'exp6_pairing_rebuild_check') && cfg.exp6_pairing_rebuild_check
    vi=numV; ti=find(update_ms>0,1,'last');
    if isempty(ti)
        o.PairingRebuildError=0;
    else
        chk=build_paired_moved_snapshot(cfg,scen0,channel_seed,x0,y0,direction, ...
            target_direction,speed_kmh(vi)/3.6,update_ms(ti)/1000);
        o.PairingRebuildError=scenario_channel_relative_error(stateBank{vi,ti},chk);
    end
else
    o.PairingRebuildError=NaN;
end

r0=evaluate_ret_ris_schemes(cfg,scen0,gamma,P_total);
if ~r0.Feasible(5)||isempty(r0.Phi{5})||~isfinite(r0.Theta(5))
    for vi=1:numV
        for ti=1:numTau
            if ~isempty(progressQueue); send(progressQueue,1); end
        end
    end
    return;
end

o.InitialFeasible=true; theta0=r0.Theta(5); phi0=r0.Phi{5}; o.InitialTheta=theta0;

gamma_stress=10^(cfg.exp6_stress_gamma_dB/10);
for vi=1:numV
    for ti=1:numTau
        scen_new=stateBank{vi,ti};
        if update_ms(ti)==0 && cfg.exp6_zero_ms_calibration
            % Definition-level calibration: no elapsed time, no passive aging.
            stale=evaluate_fixed_ret_ris_config(cfg,scen_new,theta0,phi0,gamma,P_total);
            fresh=stale; phi_fresh=phi0;
        else
            stale=evaluate_fixed_ret_ris_config(cfg,scen_new,theta0,phi0,gamma,P_total);
            % No RNG call is allowed here.  The fresh benchmark consists only
            % of deterministic warm/current-snapshot candidates and exact
            % objective acceptance, so every grid cell is reproducible.
            [phi_fresh,fresh]=optimize_fixed_ret_ris_radar(cfg,scen_new,theta0,phi0, ...
                gamma,P_total);
            if stale.feasible && (~fresh.feasible || fresh.obj_snr < stale.obj_snr*(1-1e-10))
                fresh=stale; phi_fresh=phi0; o.FreshDominanceFallback(vi,ti)=true;
            end
        end

        stress_stale=evaluate_fixed_ret_ris_config(cfg,scen_new,theta0,phi0, ...
            gamma_stress,cfg.exp6_stress_P_total);
        stress_fresh=evaluate_fixed_ret_ris_config(cfg,scen_new,theta0,phi_fresh, ...
            gamma_stress,cfg.exp6_stress_P_total);
        o.StressStaleFeasible(vi,ti)=stress_stale.feasible;
        o.StressFreshFeasible(vi,ti)=stress_fresh.feasible;
        o.PhiDistance(vi,ti)=1-abs(phi0(:)'*phi_fresh(:))/numel(phi0);
        o.PhiDistance(vi,ti)=max(0,min(1,real(o.PhiDistance(vi,ti))));
        o.StaleFeasible(vi,ti)=stale.feasible; o.FreshFeasible(vi,ti)=fresh.feasible;
        if stale.feasible; o.ObjStale(vi,ti)=stale.obj_snr; o.PcommStale(vi,ti)=stale.p_comm; end
        if fresh.feasible; o.ObjFresh(vi,ti)=fresh.obj_snr; o.PcommFresh(vi,ti)=fresh.p_comm; end
        if stale.feasible&&fresh.feasible&&stale.p_comm>0&&fresh.p_comm>0
            ppen=10*log10(stale.p_comm/fresh.p_comm);
            if abs(ppen)<=cfg.exp6_numerical_loss_tol_dB; ppen=0; end
            o.CommChange_dB(vi,ti)=ppen;
        end
        if stale.feasible&&fresh.feasible&&stale.obj_snr>0&&fresh.obj_snr>0
            raw_loss=10*log10(fresh.obj_snr/stale.obj_snr);
            if abs(raw_loss)<=cfg.exp6_numerical_loss_tol_dB; raw_loss=0; end
            o.RawLoss_dB(vi,ti)=raw_loss; o.Loss_dB(vi,ti)=raw_loss;
        end
        if ~isempty(progressQueue); send(progressQueue,1); end
    end
end
end

function scen_new=build_paired_moved_snapshot(cfg,scen0,channel_seed,x0,y0, ...
    direction,target_direction,v_mps,tau_s)
x_new=x0+direction*v_mps*tau_s;
target_x_new=cfg.target_x;
if cfg.exp6_move_sensing_targets
    target_v=cfg.exp6_target_speed_scale*v_mps;
    target_x_new=cfg.target_x+target_direction*target_v*tau_s;
end
geom_opts=struct('target_x',target_x_new,'target_y',cfg.target_y);
scen_new=generate_isac_scenario(cfg,cfg.N,cfg.K,channel_seed,x_new,y0,geom_opts);
if cfg.exp6_apply_relative_path_phase
    scen_new=apply_relative_motion_phase(scen_new,scen0,cfg.lambda);
end
end

function e=scenario_channel_relative_error(a,b)
fields={'Hd_user','Hr_user','Hd_tar','Hr_tar','G'}; e=0;
for ii=1:numel(fields)
    x=a.(fields{ii}); y=b.(fields{ii});
    e=max(e,norm(x(:)-y(:))/max(norm(x(:)),1e-15));
end
end

function sig=exp6_checkpoint_signature(cfg,speed_kmh,update_ms,gamma,P_total)
sig=sprintf('%s|%s|%s|N%d|K%d|T%d|MC%d|g%.12g|P%.12g|V%s|Tau%s|seed%d|stride%d', ...
    cfg.paper_release_version,cfg.exp6_pairing_revision,cfg.exp6_pairing_version,cfg.N,cfg.K,cfg.T,cfg.ITER_EXP6, ...
    gamma,P_total,mat2str(speed_kmh),mat2str(update_ms),cfg.seed_base,cfg.seed_stride);
end

function CP=exp6_checkpoint_struct(signature,Loss_dB,RawLoss_dB,ObjStale,ObjFresh, ...
    PcommStale,PcommFresh,CommChange_dB,PhiDistance,StaleFeasible,FreshFeasible, ...
    StressStaleFeasible,StressFreshFeasible,FreshDominanceFallback,InitialFeasible, ...
    InitialTheta,InitialX,InitialY,UserMotionDirection,TargetMotionDirection, ...
    ChannelSeed,PairingRebuildError,CompletedMC)
CP=struct(); CP.signature=signature;
CP.Loss_dB=Loss_dB; CP.RawLoss_dB=RawLoss_dB; CP.ObjStale=ObjStale; CP.ObjFresh=ObjFresh;
CP.PcommStale=PcommStale; CP.PcommFresh=PcommFresh; CP.CommChange_dB=CommChange_dB;
CP.PhiDistance=PhiDistance; CP.StaleFeasible=StaleFeasible; CP.FreshFeasible=FreshFeasible;
CP.StressStaleFeasible=StressStaleFeasible; CP.StressFreshFeasible=StressFreshFeasible;
CP.FreshDominanceFallback=FreshDominanceFallback; CP.InitialFeasible=InitialFeasible;
CP.InitialTheta=InitialTheta; CP.InitialX=InitialX; CP.InitialY=InitialY;
CP.UserMotionDirection=UserMotionDirection; CP.TargetMotionDirection=TargetMotionDirection;
CP.ChannelSeed=ChannelSeed; CP.PairingRebuildError=PairingRebuildError; CP.CompletedMC=CompletedMC;
end

function T=trend_audit_table(S,speeds,update_ms)
if isempty(fieldnames(S))
    T=table(); return;
end
rows=cell(0,15);
for v=1:numel(speeds)
    for t=1:(numel(update_ms)-1)
        rows(end+1,:)={speeds(v),update_ms(t),update_ms(t+1), ...
            S.MeanDelta(v,t),S.MeanDeltaLo(v,t),S.MeanDeltaHi(v,t), ...
            S.MedianDelta(v,t),S.NondecreasePct(v,t),S.PairCount(v,t), ...
            S.QuantileDelta(v,t),S.QuantileDeltaLo(v,t),S.QuantileDeltaHi(v,t), ...
            S.SignificantMeanDecrease(v,t),S.SignificantQuantileDecrease(v,t), ...
            (speeds(v)/3.6)*(update_ms(t+1)-update_ms(t))/1000}; %#ok<AGROW>
    end
end
T=cell2table(rows,'VariableNames',{'Speed_kmh','Tau_left_ms','Tau_right_ms', ...
    'Mean_delta_dB','Mean_delta_CI95_low_dB','Mean_delta_CI95_high_dB', ...
    'Median_delta_dB','Pairwise_nondecrease_percent','Paired_sample_count', ...
    'P90_delta_dB','P90_delta_paired_boot_CI95_low_dB','P90_delta_paired_boot_CI95_high_dB', ...
    'Significant_mean_decrease','Significant_P90_decrease','Incremental_displacement_m'});
end

function write_pairing_audit(file,cfg,ITER,speeds,update_ms,LossCount,LossP90Count, ...
    rebuildErr,fallback,TrendAudit,zeroErr)
fid=fopen(file,'w');
if fid<0; warning('Could not write pairing audit: %s',file); return; end
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'EXP.6 STRICT PAIRED MOBILITY / CHANNEL-AGING AUDIT\n');
fprintf(fid,'Release: %s\n',cfg.paper_release_version);
fprintf(fid,'Exp.6 pairing revision: %s\n',cfg.exp6_pairing_revision);
fprintf(fid,'Pairing version: %s\n',cfg.exp6_pairing_version);
fprintf(fid,'MC: %d\n',ITER);
fprintf(fid,'Speeds km/h: %s\n',mat2str(speeds));
fprintf(fid,'RIS update periods ms: %s\n',mat2str(update_ms));
fprintf(fid,'Common-random-number state bank: YES\n');
fprintf(fid,'Fresh-RIS random restarts: %d\n',logical(cfg.exp6_fresh_allow_random));
fprintf(fid,'Paired bootstrap: %d (%d replicates)\n',logical(cfg.exp6_use_paired_bootstrap),cfg.exp6_bootstrap_reps);
fprintf(fid,'Minimum paired feasible sample count: %d\n',min(LossCount(:)));
fprintf(fid,'Minimum P90 sample count: %d\n',min(LossP90Count(:)));
fprintf(fid,'Maximum state-bank rebuild relative error: %.3e\n',finite_max(rebuildErr));
fprintf(fid,'Zero-ms calibration max abs raw loss: %.3e dB\n',zeroErr);
fprintf(fid,'Fresh-dominance safeguard activations: %d\n',nnz(fallback));
if isfield(TrendAudit,'NumPointMeanDecreases')
    fprintf(fid,'Adjacent mean point-estimate decreases: %d\n',TrendAudit.NumPointMeanDecreases);
    fprintf(fid,'Adjacent mean decreases significant at paired 95%% CI: %d\n',TrendAudit.NumSignificantMeanDecreases);
    fprintf(fid,'Adjacent P90 point-estimate decreases: %d\n',TrendAudit.NumPointQuantileDecreases);
    fprintf(fid,'Adjacent P90 decreases significant at paired-bootstrap 95%% CI: %d\n',TrendAudit.NumSignificantQuantileDecreases);
end
fprintf(fid,['Interpretation: local nonmonotonicity is never clipped, smoothed, or isotonic-regressed. ', ...
    'A local reversal is treated as statistically resolved only if its paired 95%% interval lies below zero.\n']);
end

function Out=conditional_outage(Fresh,Stale)
[~,V,T]=size(Fresh); Out=nan(V,T);
for v=1:V
    for t=1:T
        fresh=Fresh(:,v,t); den=sum(fresh);
        if den>0; Out(v,t)=100*sum(fresh & ~Stale(:,v,t))/den; end
    end
end
end

function [mu,ci,n]=mean_ci_dim1(X)
[~,V,T]=size(X); mu=nan(V,T); ci=nan(V,T); n=zeros(V,T);
for v=1:V
    for t=1:T
        z=X(:,v,t); z=z(isfinite(z)); n(v,t)=numel(z);
        if isempty(z); continue; end
        mu(v,t)=mean(z);
        if numel(z)>1; ci(v,t)=student_t_critical(0.95,numel(z)-1)*std(z,0)/sqrt(numel(z)); else; ci(v,t)=0; end
    end
end
end

function q=quantile_dim1(X,p)
[~,V,T]=size(X); q=nan(V,T);
for v=1:V
    for t=1:T
        z=X(:,v,t); z=z(isfinite(z));
        if isempty(z); continue; end
        q(v,t)=empirical_quantile(z,p);
    end
end
end

function add_heatmap_values(xvals,yvals,Z,fmt)
for r=1:size(Z,1)
    for c=1:size(Z,2)
        if isfinite(Z(r,c))
            text(xvals(c),yvals(r),sprintf(fmt,Z(r,c)), ...
                'HorizontalAlignment','center','FontWeight','bold','Color','k');
        else
            text(xvals(c),yvals(r),'--','HorizontalAlignment','center','Color','k');
        end
    end
end
end

function m=finite_max(X)
z=X(isfinite(X));
if isempty(z); m=0; else; m=max(z); end
end

function m=finite_abs_max(X)
z=abs(X(isfinite(X)));
if isempty(z); m=0; else; m=max(z); end
end

function names=make_valid_names(values,prefix)
names=cell(1,numel(values));
for ii=1:numel(values)
    raw=sprintf('%s%g',prefix,values(ii));
    raw=strrep(raw,'.','p'); raw=strrep(raw,'-','m');
    names{ii}=matlab.lang.makeValidName(raw);
end
end

function save_figure_auto(fig,folder,name)
if ~exist(folder,'dir'); mkdir(folder); end
set(fig,'Color','w'); drawnow;
png=fullfile(folder,[name '.png']); figfile=fullfile(folder,[name '.fig']);
try
    exportgraphics(fig,png,'Resolution',300);
catch
    print(fig,png,'-dpng','-r300');
end
try
    savefig(fig,figfile);
catch ME
    warning('Could not save %s: %s',figfile,ME.message);
end
fprintf('Saved figure: %s\n',png);
end

function exp6_progress_bar(message)
persistent total completed t0 lastlen;
if isstruct(message)
    if isfield(message,'action')&&strcmp(message.action,'reset')
        total=max(1,double(message.total)); completed=0; t0=tic; lastlen=0;
        print_line(false);
    end
    return;
end
if isempty(total); return; end
completed=min(total,completed+double(message));
print_line(completed>=total);
if completed>=total; total=[]; completed=[]; t0=[]; lastlen=[]; end

    function print_line(done)
        elapsed=toc(t0);
        if completed>0; eta=elapsed/completed*(total-completed); else; eta=NaN; end
        L=30; filled=round(L*completed/total);
        bartext=[repmat('=',1,filled),repmat(' ',1,L-filled)];
        line=sprintf('Progress: [%s] %d/%d (%5.1f%%)  Elapsed %s  ETA %s', ...
            bartext,completed,total,100*completed/total, ...
            fmt_time(elapsed),fmt_time(eta));
        if length(line)<lastlen; line=[line repmat(' ',1,lastlen-length(line))]; end
        lastlen=length(line);
        fprintf('\r%s',line); drawnow limitrate;
        if done; fprintf('\n'); end
    end
end

function t=fmt_time(sec)
if ~isfinite(sec)||sec<0; t='--:--:--'; return; end
t=sprintf('%02d:%02d:%02d',floor(sec/3600), ...
    floor(mod(sec,3600)/60),floor(mod(sec,60)));
end
