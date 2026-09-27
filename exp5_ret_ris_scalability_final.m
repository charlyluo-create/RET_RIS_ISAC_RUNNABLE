% =========================================================================
% Experiment 5 REVISED: RET-RIS Scalability and Observed Complexity
%
% Part A: sweep RIS size N while fixing K=cfg.K.
% Part B: sweep communication-user count K while fixing N=cfg.N.
%
% Representative schemes:
%   1) Fixed RET, No RIS
%   2) RET-RIS-MinPc
%   3) RET-RIS-RadarSNR
%
% Revision principles
%   1) Never alter/post-process performance data to force a desired trend.
%   2) Use nested channels and continuation warm starts to reduce unrelated
%      Monte-Carlo variation and dimension-dependent local-optimum bias.
%   3) Reuse continuation seeds and skip preliminary Schemes 2--3 when they
%      cannot affect the reported Schemes 1/4/5. The random stream is preserved.
%   4) Extra restarts can evaluate only the two joint schemes after one full
%      initializer pass, while retaining exact objective-based acceptance.
%   5) Do not repeat the complete experiment for timing by default. The
%      separate serial timing benchmark is optional supplementary analysis;
%      Exp.7/12 provide the manuscript-level algorithm runtime comparison.
%
% Experiments 1--4 and 6 are not modified by this script.
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
ITER=cfg.ITER_EXP5;
N_range=cfg.N_range_exp5;
K_range=cfg.K_range_exp5;
gamma_dB=cfg.gamma_exp5_dB;
gamma=10^(gamma_dB/10);
P_total_dBm=cfg.P_total_exp5_dBm;
P_total=cfg.P_total_exp5;
numN=numel(N_range);
numK=numel(K_range);
num_workers=min(cfg.num_workers,ITER);
profile_env=getenv('RET_RIS_EXP5_PROFILE');
if isempty(profile_env)
    exp5_profile=cfg.exp5_execution_profile;
else
    exp5_profile=lower(strtrim(profile_env));
end
speed=resolve_exp5_speed_profile(cfg,exp5_profile);
quality_restarts=speed.quality_restarts;
use_continuation=speed.use_continuation;
run_timing=speed.run_timing;
timing_mc=speed.timing_mc;
timing_repeats=speed.timing_repeats;

% The old V4/V5.1 code repeated every N/K point in a second serial pass.
% That optional pass can dominate the total wall-clock time because the main
% Monte-Carlo simulation is parallel while the timing pass is deliberately
% serial. It is OFF by default in the final submission release. Enable it explicitly only
% when a supplementary runtime-versus-size figure is required.
force_timing_env=getenv('RET_RIS_EXP5_FORCE_TIMING');
force_timing=any(strcmpi(strtrim(force_timing_env),{'1','true','yes'}));
if force_timing
    run_timing=true;
    % Lightweight default: one complete serial sweep. Users can override.
    timing_mc=1;
    timing_repeats=1;
    tmp=str2double(getenv('RET_RIS_EXP5_TIMING_MC'));
    if isfinite(tmp) && tmp>=1; timing_mc=max(1,round(tmp)); end
    tmp=str2double(getenv('RET_RIS_EXP5_TIMING_REPEATS'));
    if isfinite(tmp) && tmp>=1; timing_repeats=max(1,round(tmp)); end
end
skip_timing_env=getenv('RET_RIS_EXP5_SKIP_TIMING');
if any(strcmpi(strtrim(skip_timing_env),{'1','true','yes'}))
    run_timing=false; % explicit skip always has highest priority
end

% Map from the unified five-scheme evaluator.
select_idx=[1 4 5];
S=numel(select_idx);
labels={'Fixed RET, No RIS','RET-RIS-MinPc','RET-RIS-RadarSNR'};
markers={'-o','-s','-p'};

%% ======================== Paths / Checks =================================
root=fileparts(mfilename('fullpath'));
if isempty(root); root=pwd; end
addpath(root);
func_path=fullfile(root,'function');
if exist(func_path,'dir'); addpath(func_path); end

required={'generate_isac_scenario','evaluate_ret_ris_schemes', ...
    'generate_channel','solve_comm_min_power','build_radar_matrix', ...
    'build_effective_user_channel','allocate_radar_power'};
for ii=1:numel(required)
    if exist(required{ii},'file')~=2
        error('Required function not found: %s.m',required{ii});
    end
end
if exist('cvx_begin','file')~=2
    error('CVX is required. Install CVX and run cvx_setup first.');
end
requested_cvx_solver=strtrim(getenv('RET_RIS_CVX_SOLVER'));
if ~isempty(requested_cvx_solver)
    try
        cvx_solver(requested_cvx_solver);
        fprintf('Requested CVX solver: %s\n',requested_cvx_solver);
    catch ME
        warning('Unable to select CVX solver %s: %s', ...
            requested_cvx_solver,ME.message);
        requested_cvx_solver='';
    end
end

stamp=datestr(now,'yyyymmdd_HHMMSS');
figure_dir=fullfile(root,sprintf('Exp5_Figures_%s_%s',lower(RUN_MODE),stamp));
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
        if ~isempty(requested_cvx_solver)
            solver_cmd=sprintf('try, cvx_solver(''%s''); catch, end', ...
                requested_cvx_solver);
            pctRunOnAll(solver_cmd);
        end
    catch ME
        warning('Parallel setup failed: %s. Falling back to serial.',ME.message);
        use_parallel=false;
    end
end

fprintf('\n=================================================================\n');
fprintf(' EXPERIMENT 5: ROBUST RET-RIS SCALABILITY / COMPLEXITY\n');
fprintf(' Mode=%s, performance MC=%d, gamma=%.1f dB, P_T=%.1f dBm\n', ...
    upper(RUN_MODE),ITER,gamma_dB,P_total_dBm);
fprintf(' N sweep: '); fprintf('%d ',N_range); fprintf('\n');
fprintf(' K sweep: '); fprintf('%d ',K_range); fprintf('\n');
fprintf(' Exp.5 execution profile=%s\n',upper(speed.name));
fprintf([' Exact-objective restarts=%d, continuation=%d, ' ...
    'selective extra restarts=%d\n'], ...
    quality_restarts,use_continuation,speed.selective_restarts);
fprintf(' Skip unused preliminary schemes on continuation points=%d\n', ...
    speed.skip_preliminary_on_continuation);
fprintf(' Reuse N-independent no-RIS baseline across N sweep=%d\n', ...
    speed.reuse_n_sweep_no_ris_baseline);
fprintf(' Continuation phase trials=%d, max AO iterations=%d\n', ...
    speed.continuation_phase_trials,speed.continuation_max_ao_iter);
if run_timing
    fprintf(' Optional timing: %d MC x %d repeats, serial after CVX warm-up\n', ...
        timing_mc,timing_repeats);
end
if SAVE_FIGURES; fprintf(' Figures: %s\n',figure_dir); end
if ~run_timing
    fprintf(' Separate timing benchmark: OFF (no second N/K solver pass)\n');
end
fprintf('=================================================================\n\n');

%% ======================== Performance Buffers ============================
ObjN=nan(ITER,numN,S);
PcommN=nan(ITER,numN,S);
FeaN=false(ITER,numN,S);
BestRestartN=nan(ITER,numN,2); % columns correspond to schemes 4 and 5

ObjK=nan(ITER,numK,S);
PcommK=nan(ITER,numK,S);
FeaK=false(ITER,numK,S);
BestRestartK=nan(ITER,numK,2);

%% ======================== Progress =======================================
progressQueue=[];
if SHOW_PROGRESS
    try
        progressQueue=parallel.pool.DataQueue;
        afterEach(progressQueue,@exp5_progress_bar);
        exp5_progress_bar(struct('action','reset','total',ITER*(numN+numK)));
    catch ME
        warning('Progress bar initialization failed: %s',ME.message);
        progressQueue=[];
    end
end

%% ======================== Performance Monte Carlo ========================
t_perf=tic;
if use_parallel
    parfor iter=1:ITER
        seed=cfg.seed_base+iter*cfg.seed_stride;
        o=run_one_mc_scalability_quality(cfg,N_range,K_range,gamma,P_total, ...
            select_idx,seed,quality_restarts,use_continuation,speed,progressQueue);
        ObjN(iter,:,:)=reshape(o.ObjN,[1 numN S]);
        PcommN(iter,:,:)=reshape(o.PcommN,[1 numN S]);
        FeaN(iter,:,:)=reshape(o.FeaN,[1 numN S]);
        BestRestartN(iter,:,:)=reshape(o.BestRestartN,[1 numN 2]);
        ObjK(iter,:,:)=reshape(o.ObjK,[1 numK S]);
        PcommK(iter,:,:)=reshape(o.PcommK,[1 numK S]);
        FeaK(iter,:,:)=reshape(o.FeaK,[1 numK S]);
        BestRestartK(iter,:,:)=reshape(o.BestRestartK,[1 numK 2]);
    end
else
    for iter=1:ITER
        seed=cfg.seed_base+iter*cfg.seed_stride;
        o=run_one_mc_scalability_quality(cfg,N_range,K_range,gamma,P_total, ...
            select_idx,seed,quality_restarts,use_continuation,speed,progressQueue);
        ObjN(iter,:,:)=reshape(o.ObjN,[1 numN S]);
        PcommN(iter,:,:)=reshape(o.PcommN,[1 numN S]);
        FeaN(iter,:,:)=reshape(o.FeaN,[1 numN S]);
        BestRestartN(iter,:,:)=reshape(o.BestRestartN,[1 numN 2]);
        ObjK(iter,:,:)=reshape(o.ObjK,[1 numK S]);
        PcommK(iter,:,:)=reshape(o.PcommK,[1 numK S]);
        FeaK(iter,:,:)=reshape(o.FeaK,[1 numK S]);
        BestRestartK(iter,:,:)=reshape(o.BestRestartK,[1 numK 2]);
    end
end
performance_runtime=toc(t_perf);

%% ======================== Robust Timing Benchmark ========================
% Timing is intentionally serial: parallel worker contention measures the
% machine scheduler rather than the algorithm. The same ascending continuation
% path and restart count used for the reported performance are timed.
TimeN_block_samples=[]; TimeN_cum_samples=[];
TimeK_block_samples=[]; TimeK_cum_samples=[];
timing_runtime=0;
if run_timing
    % Release parallel workers so the timing pass measures the algorithm on
    % one execution stream rather than resource contention across workers.
    % Parallel Computing Toolbox is optional. Release an existing pool only
    % when GCP is available; otherwise continue with the serial timing pass.
    try
        pool_for_timing=gcp('nocreate');
        if ~isempty(pool_for_timing); delete(pool_for_timing); end
    catch
        pool_for_timing=[]; %#ok<NASGU>
    end
    fprintf('\nStarting separate timing benchmark...\n');
    t_timing=tic;
    [TimeN_block_samples,TimeN_cum_samples, ...
     TimeK_block_samples,TimeK_cum_samples]=run_exp5_timing_benchmark( ...
        cfg,N_range,K_range,gamma,P_total,select_idx,timing_mc,timing_repeats, ...
        quality_restarts,use_continuation,speed);
    timing_runtime=toc(t_timing);
end

%% ======================== Performance Statistics =========================
CommonN=all(FeaN,3); CommonNRate=100*mean(CommonN,1).';
CommonK=all(FeaK,3); CommonKRate=100*mean(CommonK,1).';
[ObjN_dB,ObjN_lo,ObjN_hi,ObjN_n]=masked_mean_ci_to_db(ObjN,CommonN,1);
[PcommN_dBm,PcommN_lo,PcommN_hi,PcommN_n]= ...
    masked_mean_ci_to_db(PcommN,CommonN,1e-3);
FeaN_rate=squeeze(100*mean(FeaN,1));

[ObjK_dB,ObjK_lo,ObjK_hi,ObjK_n]=masked_mean_ci_to_db(ObjK,CommonK,1);
[PcommK_dBm,PcommK_lo,PcommK_hi,PcommK_n]= ...
    masked_mean_ci_to_db(PcommK,CommonK,1e-3);
FeaK_rate=squeeze(100*mean(FeaK,1));

ObjNRef=repmat(ObjN(:,1,:),[1 numN 1]);
MaskNPair=FeaN & repmat(FeaN(:,1,:),[1 numN 1]);
[ObjNGain_dB,ObjNGain_lo,ObjNGain_hi,ObjNGain_n]= ...
    paired_db_statistics(ObjN,ObjNRef,MaskNPair);

% Median/IQR timing summaries. Cumulative time is the complete executable
% path; block time is saved as a diagnostic for algorithmic attribution.
[TimeN_med,TimeN_q1,TimeN_q3]=sample_summary(TimeN_cum_samples,numN,S);
[TimeN_block_med,TimeN_block_q1,TimeN_block_q3]= ...
    sample_summary(TimeN_block_samples,numN,S);
[TimeK_med,TimeK_q1,TimeK_q3]=sample_summary(TimeK_cum_samples,numK,S);
[TimeK_block_med,TimeK_block_q1,TimeK_block_q3]= ...
    sample_summary(TimeK_block_samples,numK,S);

%% ======================== Complexity Proxies =============================
% These are normalized arithmetic-work proxies, not measured time and not
% formal big-O equalities. They expose the N-linear phase-gradient workload
% and the stronger K dependence of multi-user terms.
[WorkN,WorkK]=complexity_proxies(cfg,N_range,K_range);

if SAVE_MAIN_FIGURE
%% ======================== Main Paper Composite ===========================
% The manuscript uses one 2x2 scalability figure. Detailed runtime,
% feasibility and workload plots are retained as supplementary outputs.
fig0=figure('Position',[30 30 1450 900],'Color','w');
tl0=tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

nexttile;
for s=1:S
    errorbar(N_range,ObjN_dB(:,s),ObjN_dB(:,s)-ObjN_lo(:,s), ...
        ObjN_hi(:,s)-ObjN_dB(:,s),markers{s},'LineWidth',1.5,'MarkerSize',7); hold on;
end
grid on; xlabel('RIS elements N'); ylabel('Weighted radar SNR (dB)');
title('(a) Sensing scalability with RIS size'); legend(labels,'Location','best');

nexttile;
for s=1:S
    errorbar(N_range,PcommN_dBm(:,s),PcommN_dBm(:,s)-PcommN_lo(:,s), ...
        PcommN_hi(:,s)-PcommN_dBm(:,s),markers{s},'LineWidth',1.5,'MarkerSize',7); hold on;
end
grid on; xlabel('RIS elements N'); ylabel('Communication power (dBm)');
title('(b) Communication cost with RIS size');

nexttile;
for s=1:S
    errorbar(K_range,ObjK_dB(:,s),ObjK_dB(:,s)-ObjK_lo(:,s), ...
        ObjK_hi(:,s)-ObjK_dB(:,s),markers{s},'LineWidth',1.5,'MarkerSize',7); hold on;
end
grid on; xlabel('Communication users K'); ylabel('Weighted radar SNR (dB)');
title('(c) Sensing performance with network load');

nexttile;
for s=1:S
    errorbar(K_range,PcommK_dBm(:,s),PcommK_dBm(:,s)-PcommK_lo(:,s), ...
        PcommK_hi(:,s)-PcommK_dBm(:,s),markers{s},'LineWidth',1.5,'MarkerSize',7); hold on;
end
grid on; xlabel('Communication users K'); ylabel('Communication power (dBm)');
title('(d) Communication cost with network load');

sgtitle(tl0,sprintf('RET-RIS Scalability (QoS = %.1f dB, P_T = %.1f dBm)', ...
    gamma_dB,P_total_dBm));
if SAVE_MAIN_FIGURE; save_figure_auto(fig0,figure_dir,'00_Main_Scalability_Composite'); end

end

if SAVE_SUPP_FIGURES
%% ======================== Figure 1: Radar SNR vs N ========================
fig1=figure('Position',[60 60 950 620],'Color','w');
for s=1:S
    errorbar(N_range,ObjN_dB(:,s),ObjN_dB(:,s)-ObjN_lo(:,s), ...
        ObjN_hi(:,s)-ObjN_dB(:,s),markers{s}, ...
        'LineWidth',1.6,'MarkerSize',8); hold on;
end
xlabel('Number of RIS elements N');
ylabel('Weighted radar SNR (dB)');
legend(labels,'Location','northwest'); grid on;
title(sprintf('Sensing Scalability with RIS Size (K=%d)',cfg.K));
if SAVE_SUPP_FIGURES; save_figure_auto(fig1,figure_dir,'S1_Radar_SNR_vs_RIS_Elements'); end

%% ======================== Figure 2: Pcomm vs N ============================
fig2=figure('Position',[80 80 950 620],'Color','w');
for s=1:S
    errorbar(N_range,PcommN_dBm(:,s),PcommN_dBm(:,s)-PcommN_lo(:,s), ...
        PcommN_hi(:,s)-PcommN_dBm(:,s),markers{s}, ...
        'LineWidth',1.6,'MarkerSize',8); hold on;
end
xlabel('Number of RIS elements N');
ylabel('Minimum communication power P_c^{min} (dBm)');
legend(labels,'Location','northeast'); grid on;
title(sprintf('Communication Scalability with RIS Size (K=%d)',cfg.K));
if SAVE_SUPP_FIGURES; save_figure_auto(fig2,figure_dir,'S2_Communication_Power_vs_RIS_Elements'); end

%% ======================== Figure 3: Runtime vs N ==========================
if run_timing
    fig3=figure('Position',[100 100 950 620],'Color','w');
    for s=1:S
        errorbar(N_range,TimeN_med(:,s),TimeN_med(:,s)-TimeN_q1(:,s), ...
            TimeN_q3(:,s)-TimeN_med(:,s),markers{s}, ...
            'LineWidth',1.6,'MarkerSize',8); hold on;
    end
    xlabel('Number of RIS elements N');
    ylabel('Median cumulative runtime per realization (s)');
    legend(labels,'Location','northwest'); grid on;
    title('Observed Computational Scaling after Warm-Up');
    if SAVE_SUPP_FIGURES
        save_figure_auto(fig3,figure_dir,'S3_Runtime_vs_RIS_Elements');
    end
end

%% ======================== Figure 4: Radar SNR vs K ========================
fig4=figure('Position',[120 120 950 620],'Color','w');
for s=1:S
    errorbar(K_range,ObjK_dB(:,s),ObjK_dB(:,s)-ObjK_lo(:,s), ...
        ObjK_hi(:,s)-ObjK_dB(:,s),markers{s}, ...
        'LineWidth',1.6,'MarkerSize',8); hold on;
end
xlabel('Number of communication users K');
ylabel('Weighted radar SNR (dB)');
legend(labels,'Location','southwest'); grid on;
title(sprintf('Sensing Performance under Increasing IoT Load (N=%d)',cfg.N));
if SAVE_SUPP_FIGURES; save_figure_auto(fig4,figure_dir,'S4_Radar_SNR_vs_User_Count'); end

%% ======================== Figure 5: Feasibility vs K ======================
fig5=figure('Position',[140 140 950 620],'Color','w');
for s=1:S
    plot(K_range,FeaK_rate(:,s),markers{s},'LineWidth',1.6,'MarkerSize',8); hold on;
end
xlabel('Number of communication users K');
ylabel('Feasibility rate (%)'); ylim([0 105]);
legend(labels,'Location','southwest'); grid on;
title(sprintf('Network-Load Feasibility at P_T=%.1f dBm',P_total_dBm));
if SAVE_SUPP_FIGURES; save_figure_auto(fig5,figure_dir,'S5_Feasibility_vs_User_Count'); end

%% ======================== Figure 6: Runtime vs K ==========================
if run_timing
    fig6=figure('Position',[160 160 950 620],'Color','w');
    for s=1:S
        errorbar(K_range,TimeK_med(:,s),TimeK_med(:,s)-TimeK_q1(:,s), ...
            TimeK_q3(:,s)-TimeK_med(:,s),markers{s}, ...
            'LineWidth',1.6,'MarkerSize',8); hold on;
    end
    xlabel('Number of communication users K');
    ylabel('Median cumulative runtime per realization (s)');
    legend(labels,'Location','northwest'); grid on;
    title('Observed Computational Scaling with Network Load');
    if SAVE_SUPP_FIGURES
        save_figure_auto(fig6,figure_dir,'S6_Runtime_vs_User_Count');
    end
end


%% ======================== Figure 7: Paired RIS-Size Gain =================
fig7=figure('Position',[180 180 980 630],'Color','w'); hold on;
for ss=1:S
    errorbar(N_range,ObjNGain_dB(:,ss), ...
        ObjNGain_dB(:,ss)-ObjNGain_lo(:,ss), ...
        ObjNGain_hi(:,ss)-ObjNGain_dB(:,ss),markers{ss}, ...
        'LineWidth',1.6,'MarkerSize',8);
end
yline(0,'--k','HandleVisibility','off'); grid on;
xlabel('Number of RIS elements N');
ylabel(sprintf('Paired radar-SNR gain over N=%d (dB)',N_range(1)));
legend(labels,'Location','best');
title('Nested-Channel Paired RIS-Size Gain');
if SAVE_SUPP_FIGURES; save_figure_auto(fig7,figure_dir,'S7_Paired_RIS_Size_Gain'); end

%% ======================== Figure 8: Complexity Proxy ======================
fig8=figure('Position',[200 200 1100 470],'Color','w');
subplot(1,2,1);
plot(N_range,WorkN(:,1),'-s','LineWidth',1.6,'MarkerSize',8); hold on;
plot(N_range,WorkN(:,2),'-p','LineWidth',1.6,'MarkerSize',8);
xlabel('Number of RIS elements N'); ylabel('Normalized arithmetic workload');
legend({'RET-RIS-MinPc','RET-RIS-RadarSNR'},'Location','northwest');
grid on; title('Phase-Optimization Workload versus N');
subplot(1,2,2);
plot(K_range,WorkK(:,1),'-s','LineWidth',1.6,'MarkerSize',8); hold on;
plot(K_range,WorkK(:,2),'-p','LineWidth',1.6,'MarkerSize',8);
xlabel('Number of users K'); ylabel('Normalized arithmetic workload');
legend({'RET-RIS-MinPc','RET-RIS-RadarSNR'},'Location','northwest');
grid on; title('Multiuser Workload versus K');
if SAVE_SUPP_FIGURES; save_figure_auto(fig8,figure_dir,'S8_Normalized_Complexity_Proxy'); end

end

%% ======================== Console / Save ==================================
fprintf('\nAverage radar SNR versus N (dB):\n');
disp(array2table([N_range(:),ObjN_dB], ...
    'VariableNames',{'N','Fixed_NoRIS','RET_RIS_MinPc','RET_RIS_RadarSNR'}));
fprintf('\nAverage communication power versus N (dBm):\n');
disp(array2table([N_range(:),PcommN_dBm], ...
    'VariableNames',{'N','Fixed_NoRIS','RET_RIS_MinPc','RET_RIS_RadarSNR'}));
fprintf('\nAverage radar SNR versus K (dB):\n');
disp(array2table([K_range(:),ObjK_dB], ...
    'VariableNames',{'K','Fixed_NoRIS','RET_RIS_MinPc','RET_RIS_RadarSNR'}));
fprintf('\nFeasibility versus K (%%):\n');
disp(array2table([K_range(:),FeaK_rate], ...
    'VariableNames',{'K','Fixed_NoRIS','RET_RIS_MinPc','RET_RIS_RadarSNR'}));
if run_timing
    fprintf('\nMedian cumulative runtime versus N (s):\n');
    disp(array2table([N_range(:),TimeN_med], ...
        'VariableNames',{'N','Fixed_NoRIS','RET_RIS_MinPc','RET_RIS_RadarSNR'}));
    fprintf('\nMedian cumulative runtime versus K (s):\n');
    disp(array2table([K_range(:),TimeK_med], ...
        'VariableNames',{'K','Fixed_NoRIS','RET_RIS_MinPc','RET_RIS_RadarSNR'}));
end

results_dir=fullfile(root,'results');
if ~exist(results_dir,'dir'); mkdir(results_dir); end
save(fullfile(results_dir,sprintf('exp5_scalability_revised_%s.mat',lower(RUN_MODE))), ...
    'cfg','N_range','K_range','gamma_dB','gamma','P_total','P_total_dBm', ...
    'ObjN','PcommN','FeaN','BestRestartN','CommonN','CommonNRate', ...
    'ObjNGain_dB','ObjNGain_lo','ObjNGain_hi','ObjNGain_n', ...
    'ObjK','PcommK','FeaK','BestRestartK','CommonK','CommonKRate', ...
    'ObjN_dB','PcommN_dBm','FeaN_rate','ObjK_dB','PcommK_dBm','FeaK_rate', ...
    'TimeN_block_samples','TimeN_cum_samples','TimeK_block_samples','TimeK_cum_samples', ...
    'TimeN_med','TimeN_q1','TimeN_q3','TimeN_block_med','TimeN_block_q1','TimeN_block_q3', ...
    'TimeK_med','TimeK_q1','TimeK_q3','TimeK_block_med','TimeK_block_q1','TimeK_block_q3', ...
    'WorkN','WorkK','performance_runtime','timing_runtime', ...
    'speed','exp5_profile','requested_cvx_solver','run_timing','force_timing', ...
    'timing_mc','timing_repeats');

total_time=toc(t_total);
fprintf('\n=================================================================\n');
fprintf(' EXPERIMENT 5 REVISED COMPLETED\n');
fprintf(' Performance MC: %.1f s (%.2f min / %.2f h)\n', ...
    performance_runtime,performance_runtime/60,performance_runtime/3600);
fprintf(' Timing benchmark: %.1f s (%.2f min / %.2f h)\n', ...
    timing_runtime,timing_runtime/60,timing_runtime/3600);
fprintf(' Total runtime:     %.1f s (%.2f min / %.2f h)\n', ...
    total_time,total_time/60,total_time/3600);
if SAVE_FIGURES; fprintf(' Figure folder:     %s\n',figure_dir); end
fprintf('=================================================================\n');

%% ========================================================================
%% Local functions
%% ========================================================================
function o=run_one_mc_scalability_quality(cfg,N_range,K_range,gamma,P_total, ...
    select_idx,seed,restarts,use_continuation,speed,progressQueue)

S=numel(select_idx); numN=numel(N_range); numK=numel(K_range);
o.ObjN=nan(numN,S); o.PcommN=nan(numN,S); o.FeaN=false(numN,S);
o.BestRestartN=nan(numN,2);
o.ObjK=nan(numK,S); o.PcommK=nan(numK,S); o.FeaK=false(numK,S);
o.BestRestartK=nan(numK,2);

% -------------------- Paired RIS-size sweep -------------------------------
rng(seed,'twister');
x=cfg.road_x_min+(cfg.road_x_max-cfg.road_x_min)*rand(1,cfg.K);
y=-cfg.road_y_half+2*cfg.road_y_half*rand(1,cfg.K);
scenMaxN=generate_isac_scenario(cfg,max(N_range),cfg.K,seed+100,x,y,struct());

prevPhiMin=[]; prevPhiRadar=[]; prevThetaMin=[]; prevThetaRadar=[];
baselineN=[];
for nIdx=1:numN
    scen=subset_ris_scenario(scenMaxN,N_range(nIdx));
    opts=struct();
    if nIdx>1 && speed.reuse_n_sweep_no_ris_baseline
        opts.skip_scheme1=true;
    end
    if use_continuation
        opts.phi_seed_minpc=extend_phase_seed(prevPhiMin,scen.N);
        opts.phi_seed_radar=extend_phase_seed(prevPhiRadar,scen.N);
        opts.theta_seed_minpc_deg=prevThetaMin;
        opts.theta_seed_radar_deg=prevThetaRadar;
        opts.num_phase_trials_override=speed.continuation_phase_trials;
        opts.max_ao_iter_override=speed.continuation_max_ao_iter;
        if speed.skip_preliminary_on_continuation && ...
                ~isempty(opts.phi_seed_minpc) && ...
                ~isempty(opts.theta_seed_minpc_deg)
            opts.skip_preliminary_schemes=true;
            opts.preserve_rng_when_skipping=true;
        end
    end
    [r,bestRestart]=evaluate_quality_controlled(cfg,scen,gamma,P_total, ...
        restarts,seed+10000+nIdx*503,opts,speed);
    if nIdx==1
        baselineN=r;
    elseif speed.reuse_n_sweep_no_ris_baseline && ~isempty(baselineN)
        r=copy_scheme(r,baselineN,1);
    end
    o.ObjN(nIdx,:)=r.Obj(select_idx);
    o.PcommN(nIdx,:)=r.Pcomm(select_idx);
    o.FeaN(nIdx,:)=r.Feasible(select_idx);
    o.BestRestartN(nIdx,:)=bestRestart;
    if use_continuation
        if r.Feasible(4); prevPhiMin=r.Phi{4}; prevThetaMin=r.Theta(4); end
        if r.Feasible(5); prevPhiRadar=r.Phi{5}; prevThetaRadar=r.Theta(5); end
    end
    if ~isempty(progressQueue); send(progressQueue,1); end
end

% -------------------- Nested user-load sweep ------------------------------
Kmax=max(K_range);
rng(seed+7000,'twister');
xmax=cfg.road_x_min+(cfg.road_x_max-cfg.road_x_min)*rand(1,Kmax);
ymax=-cfg.road_y_half+2*cfg.road_y_half*rand(1,Kmax);
scenMax=generate_isac_scenario(cfg,cfg.N,Kmax,seed+7100,xmax,ymax,struct());

prevPhiMin=[]; prevPhiRadar=[]; prevThetaMin=[]; prevThetaRadar=[];
for kIdx=1:numK
    scen=subset_scenario(scenMax,K_range(kIdx));
    opts=struct();
    if use_continuation
        opts.phi_seed_minpc=prevPhiMin;
        opts.phi_seed_radar=prevPhiRadar;
        opts.theta_seed_minpc_deg=prevThetaMin;
        opts.theta_seed_radar_deg=prevThetaRadar;
        opts.num_phase_trials_override=speed.continuation_phase_trials;
        opts.max_ao_iter_override=speed.continuation_max_ao_iter;
        if speed.skip_preliminary_on_continuation && ...
                ~isempty(opts.phi_seed_minpc) && ...
                ~isempty(opts.theta_seed_minpc_deg)
            opts.skip_preliminary_schemes=true;
            opts.preserve_rng_when_skipping=true;
        end
    end
    [r,bestRestart]=evaluate_quality_controlled(cfg,scen,gamma,P_total, ...
        restarts,seed+20000+kIdx*907,opts,speed);
    o.ObjK(kIdx,:)=r.Obj(select_idx);
    o.PcommK(kIdx,:)=r.Pcomm(select_idx);
    o.FeaK(kIdx,:)=r.Feasible(select_idx);
    o.BestRestartK(kIdx,:)=bestRestart;
    if use_continuation
        if r.Feasible(4); prevPhiMin=r.Phi{4}; prevThetaMin=r.Theta(4); end
        if r.Feasible(5); prevPhiRadar=r.Phi{5}; prevThetaRadar=r.Theta(5); end
    end
    if ~isempty(progressQueue); send(progressQueue,1); end
end
end

function [best,bestRestart]=evaluate_quality_controlled(cfg,scen,gamma,P_total, ...
    restarts,seed,opts,speed)
% Scheme 4 is selected by exact minimum Pcomm; scheme 5 by exact maximum
% radar SNR. No monotonic envelope or synthetic correction is applied.
best=[]; bestRestart=[NaN NaN];
runtimeBlockTotal=zeros(1,5); runtimeCumulativeTotal=zeros(1,5);
for rr=1:restarts
    rng(seed+rr*104729,'twister');
    opts_rr=opts;

    % The first restart computes the complete five-scheme initializer path.
    % Extra restarts can reuse the current point's exact feasible joint
    % solutions and run only Schemes 4--5. This removes repeated evaluation
    % of the deterministic fixed baseline and exhaustive RET-only scan.
    can_selective=rr>1 && speed.selective_restarts && ~isempty(best) && ...
        best.Feasible(4) && best.Feasible(5);
    if can_selective
        opts_rr.phi_seed_minpc=best.Phi{4};
        opts_rr.theta_seed_minpc_deg=best.Theta(4);
        opts_rr.phi_seed_radar=best.Phi{5};
        opts_rr.theta_seed_radar_deg=best.Theta(5);
        opts_rr.skip_scheme1=true;
        opts_rr.skip_preliminary_schemes=true;
        opts_rr.preserve_rng_when_skipping=true;
    end

    cand=evaluate_ret_ris_schemes(cfg,scen,gamma,P_total,opts_rr);
    rb=cand.RuntimeBlock; rb(~isfinite(rb))=0;
    rc=cand.RuntimeCumulative; rc(~isfinite(rc))=0;
    runtimeBlockTotal=runtimeBlockTotal+rb;
    runtimeCumulativeTotal=runtimeCumulativeTotal+rc;
    if isempty(best)
        best=cand;
        if best.Feasible(4); bestRestart(1)=rr; end
        if best.Feasible(5); bestRestart(2)=rr; end
    else
        if cand.Feasible(4) && (~best.Feasible(4) || cand.Pcomm(4)<best.Pcomm(4))
            best=copy_scheme(best,cand,4); bestRestart(1)=rr;
        end
        if cand.Feasible(5) && (~best.Feasible(5) || cand.Obj(5)>best.Obj(5))
            best=copy_scheme(best,cand,5); bestRestart(2)=rr;
        end
    end
end
if ~isempty(best)
    best.RuntimeBlock=runtimeBlockTotal;
    best.RuntimeCumulative=runtimeCumulativeTotal;
end
end

function a=copy_scheme(a,b,s)
fields={'Pcomm','Pradar','Obj','Feasible','Theta','G_eff','Jcomm','Jradar', ...
    'RuntimeBlock','RuntimeCumulative'};
for ii=1:numel(fields)
    f=fields{ii}; a.(f)(s)=b.(f)(s);
end
a.Phi{s}=b.Phi{s};
end

function phiNew=extend_phase_seed(phiOld,Nnew)
% ULA-aware continuation: preserve existing entries and extrapolate their
% average phase increment. It is only an initializer; exact objectives decide
% whether the candidate is accepted.
if isempty(phiOld); phiNew=[]; return; end
phiOld=phiOld(:); phiOld=phiOld./max(abs(phiOld),eps);
Nold=numel(phiOld);
if Nnew<=Nold
    phiNew=phiOld(1:Nnew); return;
end
if Nold>=2
    inc=angle(phiOld(2:end).*conj(phiOld(1:end-1)));
    delta=angle(mean(exp(1i*inc)));
else
    delta=0;
end
extra=(1:(Nnew-Nold)).';
phiNew=[phiOld; phiOld(end).*exp(1i*delta*extra)];
phiNew=phiNew./max(abs(phiNew),eps);
end

function [TimeNBlock,TimeNCum,TimeKBlock,TimeKCum]= ...
    run_exp5_timing_benchmark(cfg,N_range,K_range,gamma,P_total, ...
    select_idx,timing_mc,timing_repeats,restarts,use_continuation,speed)
% Serial timing with the same restarts and continuation as the performance
% pass. Each timing sample is one complete ascending N sweep and one complete
% ascending K sweep. Reported times therefore correspond to the performance
% actually shown, rather than a cheaper single-start surrogate.

numN=numel(N_range); numK=numel(K_range); S=numel(select_idx);
numSamples=timing_mc*timing_repeats;
TimeNBlock=nan(numSamples,numN,S); TimeNCum=nan(numSamples,numN,S);
TimeKBlock=nan(numSamples,numK,S); TimeKCum=nan(numSamples,numK,S);

fprintf('  Timing warm-up...\n');
for ww=1:max(1,cfg.exp5_timing_warmup_runs)
    seedWarm=cfg.seed_base+900000+ww;
    scenWarm=generate_isac_scenario(cfg,min(N_range),min(K_range),seedWarm,[],[],struct());
    cfgWarm=cfg; cfgWarm.max_ao_iter=min(2,cfg.max_ao_iter);
    cfgWarm.num_phase_trials=min(4,cfg.num_phase_trials);
    evaluate_quality_controlled(cfgWarm,scenWarm,gamma,P_total,1,seedWarm+1,struct(),speed);
    try; cvx_clear; catch; end
end

sample=0;
for mc=1:timing_mc
    baseSeed=cfg.seed_base+800000+mc*cfg.seed_stride;
    for rep=1:timing_repeats
        sample=sample+1;
        seed=baseSeed+rep*3001;

        rng(seed,'twister');
        x=cfg.road_x_min+(cfg.road_x_max-cfg.road_x_min)*rand(1,cfg.K);
        y=-cfg.road_y_half+2*cfg.road_y_half*rand(1,cfg.K);
        scenMaxN=generate_isac_scenario(cfg,max(N_range),cfg.K,seed+100,x,y,struct());
        prevPhiMin=[]; prevPhiRadar=[]; prevThetaMin=[]; prevThetaRadar=[];
        baselineN=[];
        for nIdx=1:numN
            scen=subset_ris_scenario(scenMaxN,N_range(nIdx));
            opts=struct();
            if nIdx>1 && speed.reuse_n_sweep_no_ris_baseline
                opts.skip_scheme1=true;
            end
            if use_continuation
                opts.phi_seed_minpc=extend_phase_seed(prevPhiMin,scen.N);
                opts.phi_seed_radar=extend_phase_seed(prevPhiRadar,scen.N);
                opts.theta_seed_minpc_deg=prevThetaMin;
                opts.theta_seed_radar_deg=prevThetaRadar;
                opts.num_phase_trials_override=speed.continuation_phase_trials;
                opts.max_ao_iter_override=speed.continuation_max_ao_iter;
                if speed.skip_preliminary_on_continuation && ...
                        ~isempty(opts.phi_seed_minpc) && ...
                        ~isempty(opts.theta_seed_minpc_deg)
                    opts.skip_preliminary_schemes=true;
                    opts.preserve_rng_when_skipping=true;
                end
            end
            [r,~]=evaluate_quality_controlled(cfg,scen,gamma,P_total, ...
                restarts,seed+10000+nIdx*503,opts,speed);
            if nIdx==1
                baselineN=r;
            elseif speed.reuse_n_sweep_no_ris_baseline && ~isempty(baselineN)
                r=copy_scheme(r,baselineN,1);
            end
            TimeNBlock(sample,nIdx,:)=reshape(r.RuntimeBlock(select_idx),[1 1 S]);
            TimeNCum(sample,nIdx,:)=reshape(r.RuntimeCumulative(select_idx),[1 1 S]);
            if use_continuation
                if r.Feasible(4); prevPhiMin=r.Phi{4}; prevThetaMin=r.Theta(4); end
                if r.Feasible(5); prevPhiRadar=r.Phi{5}; prevThetaRadar=r.Theta(5); end
            end
        end

        Kmax=max(K_range);
        rng(seed+7000,'twister');
        xmax=cfg.road_x_min+(cfg.road_x_max-cfg.road_x_min)*rand(1,Kmax);
        ymax=-cfg.road_y_half+2*cfg.road_y_half*rand(1,Kmax);
        scenMaxK=generate_isac_scenario(cfg,cfg.N,Kmax,seed+7100,xmax,ymax,struct());
        prevPhiMin=[]; prevPhiRadar=[]; prevThetaMin=[]; prevThetaRadar=[];
        for kIdx=1:numK
            scen=subset_scenario(scenMaxK,K_range(kIdx));
            opts=struct();
            if use_continuation
                opts.phi_seed_minpc=prevPhiMin;
                opts.phi_seed_radar=prevPhiRadar;
                opts.theta_seed_minpc_deg=prevThetaMin;
                opts.theta_seed_radar_deg=prevThetaRadar;
                opts.num_phase_trials_override=speed.continuation_phase_trials;
                opts.max_ao_iter_override=speed.continuation_max_ao_iter;
                if speed.skip_preliminary_on_continuation && ...
                        ~isempty(opts.phi_seed_minpc) && ...
                        ~isempty(opts.theta_seed_minpc_deg)
                    opts.skip_preliminary_schemes=true;
                    opts.preserve_rng_when_skipping=true;
                end
            end
            [r,~]=evaluate_quality_controlled(cfg,scen,gamma,P_total, ...
                restarts,seed+20000+kIdx*907,opts,speed);
            TimeKBlock(sample,kIdx,:)=reshape(r.RuntimeBlock(select_idx),[1 1 S]);
            TimeKCum(sample,kIdx,:)=reshape(r.RuntimeCumulative(select_idx),[1 1 S]);
            if use_continuation
                if r.Feasible(4); prevPhiMin=r.Phi{4}; prevThetaMin=r.Theta(4); end
                if r.Feasible(5); prevPhiRadar=r.Phi{5}; prevThetaRadar=r.Theta(5); end
            end
        end
        fprintf('  Timing sample %d/%d completed.\n',sample,numSamples);
    end
end
end

function scen=subset_ris_scenario(scenMax,N)
scen=scenMax;
if N>scenMax.N; error('Requested N exceeds generated maximum.'); end
scen.N=N;
scen.Hr_user=scenMax.Hr_user(1:N,:);
scen.Hr_tar=scenMax.Hr_tar(1:N,:);
scen.G=scenMax.G(1:N,:);
end

function scen=subset_scenario(scenMax,K)
scen=scenMax;
if K>scenMax.K; error('Requested K exceeds generated maximum.'); end
scen.K=K;
scen.x_u=scenMax.x_u(1:K);
scen.y_u=scenMax.y_u(1:K);
scen.Hd_user=scenMax.Hd_user(:,1:K);
scen.Hr_user=scenMax.Hr_user(:,1:K);
scen.gain_u=scenMax.gain_u(:,1:K);
scen.elev_users=scenMax.elev_users(1:K);
end

function [mu_db,lo_db,hi_db]=curve_linear_ci_to_db(X,reference)
[~,P,S]=size(X);
mu_db=nan(P,S); lo_db=mu_db; hi_db=mu_db;
for p=1:P
    for s=1:S
        v=X(:,p,s); v=v(isfinite(v));
        if isempty(v); continue; end
        mu=mean(v);
        if numel(v)>1; half=student_t_critical(0.95,numel(v)-1)*std(v,0)/sqrt(numel(v)); else; half=0; end
        lo=max(mu-half,eps); hi=max(mu+half,eps);
        mu_db(p,s)=10*log10(max(mu,eps)/reference);
        lo_db(p,s)=10*log10(lo/reference);
        hi_db(p,s)=10*log10(hi/reference);
    end
end
end

function [med,q1,q3]=sample_summary(X,P,S)
med=nan(P,S); q1=med; q3=med;
if isempty(X); return; end
for p=1:P
    for s=1:S
        v=X(:,p,s); v=v(isfinite(v));
        if isempty(v); continue; end
        v=sort(v(:));
        med(p,s)=quantile_manual(v,0.50);
        q1(p,s)=quantile_manual(v,0.25);
        q3(p,s)=quantile_manual(v,0.75);
    end
end
end

function q=quantile_manual(v,p)
n=numel(v);
if n==1; q=v(1); return; end
x=1+(n-1)*p;
lo=floor(x); hi=ceil(x);
if lo==hi; q=v(lo); else; q=v(lo)+(x-lo)*(v(hi)-v(lo)); end
end

function [WorkN,WorkK]=complexity_proxies(cfg,N_range,K_range)
M=cfg.M; T=cfg.T; L=cfg.max_ao_iter; R=100;
K0=cfg.K; N0=cfg.N;
workN_min=L*R*(N_range(:)*M*K0^2);
workN_rad=L*R*(N_range(:)*M*(K0^2+T));
WorkN=[workN_min/workN_min(1),workN_rad/workN_rad(1)];
workK_min=L*R*(N0*M*(K_range(:).^2));
workK_rad=L*R*(N0*M*(K_range(:).^2+T));
WorkK=[workK_min/workK_min(1),workK_rad/workK_rad(1)];
end

function speed=resolve_exp5_speed_profile(cfg,name)
% Three profiles are provided:
%   reproducible : original v5.1 execution path and timing repetitions.
%   balanced     : same PAPER Monte Carlo/restart count, but extra restarts
%                  reuse exact joint solutions and timing uses 3 full sweeps.
%   fast         : debugging/preview profile with one restart and reduced
%                  continuation search effort.
name=lower(strtrim(char(name)));
speed=struct();
switch name
    case {'reproducible','original','exact'}
        speed.name='reproducible';
        speed.quality_restarts=max(1,cfg.exp5_quality_restarts);
        speed.use_continuation=logical(cfg.exp5_use_continuation);
        speed.selective_restarts=false;
        speed.skip_preliminary_on_continuation=true;
        speed.reuse_n_sweep_no_ris_baseline=true;
        speed.continuation_phase_trials=cfg.num_phase_trials;
        speed.continuation_max_ao_iter=cfg.max_ao_iter;
        speed.run_timing=logical(cfg.exp5_run_timing_benchmark);
        if strcmpi(cfg.RUN_MODE,'paper')
            speed.timing_mc=5; speed.timing_repeats=3;
        else
            speed.timing_mc=max(1,cfg.exp5_timing_mc);
            speed.timing_repeats=max(1,cfg.exp5_timing_repeats);
        end
    case {'balanced','paper'}
        speed.name='balanced';
        speed.quality_restarts=max(1,cfg.exp5_quality_restarts);
        speed.use_continuation=logical(cfg.exp5_use_continuation);
        speed.selective_restarts=true;
        speed.skip_preliminary_on_continuation=true;
        speed.reuse_n_sweep_no_ris_baseline=true;
        speed.continuation_phase_trials=cfg.num_phase_trials;
        speed.continuation_max_ao_iter=cfg.max_ao_iter;
        speed.run_timing=logical(cfg.exp5_run_timing_benchmark);
        speed.timing_mc=max(1,cfg.exp5_timing_mc);
        speed.timing_repeats=max(1,cfg.exp5_timing_repeats);
    case {'fast','preview'}
        speed.name='fast';
        speed.quality_restarts=1;
        speed.use_continuation=true;
        speed.selective_restarts=true;
        speed.skip_preliminary_on_continuation=true;
        speed.reuse_n_sweep_no_ris_baseline=true;
        speed.continuation_phase_trials=max(4,min(6,cfg.num_phase_trials));
        speed.continuation_max_ao_iter=max(3,min(6,cfg.max_ao_iter));
        speed.run_timing=logical(cfg.exp5_run_timing_benchmark);
        speed.timing_mc=1; speed.timing_repeats=1;
    otherwise
        error(['Unknown RET_RIS_EXP5_PROFILE=%s. Use reproducible, ' ...
            'balanced, or fast.'],name);
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

function exp5_progress_bar(message)
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
