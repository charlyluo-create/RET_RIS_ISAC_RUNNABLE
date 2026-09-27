% =========================================================================
% Experiment 9 (): Imperfect-CSI Sensitivity on a Target-NMSE Grid
%
% Hardware variables and active beams are designed from estimated channels,
% then evaluated WITHOUT re-optimization on the true channels.  The error
% model is
%   Hhat = sqrt(1-epsilon^2) H + epsilon E,
% where E is independently generated and normalized to the energy of each
% link.  The manuscript axis is the target error-power ratio epsilon^2 in
% dB ("Perfect", -30, -20, -15, -10, -8, -5 dB in PAPER mode), rather than the
% amplitude coefficient epsilon.  Actual measured link-wise NMSE is saved.
%
% Primary sensing statistics use a FIXED perfect-CSI cohort.  If an estimated
% design violates the actual communication QoS, its sensing utility is set to
% zero in the outage-penalized metric.  This avoids survivor bias at high CSI
% error.  Conditional-on-success results are retained only for diagnosis.
% =========================================================================
clear; close all; clc;
t_total=tic; SHOW_PROGRESS=true;
env_mode=getenv('RET_RIS_RUN_MODE');
if isempty(env_mode); RUN_MODE='paper'; else; RUN_MODE=lower(env_mode); end
cfg=unified_ret_ris_isac_config(RUN_MODE);
FIG_POLICY=get_figure_policy(RUN_MODE);
SAVE_MAIN_FIGURE=FIG_POLICY.save_main;
SAVE_SUPP_FIGURES=FIG_POLICY.save_supplementary; %#ok<NASGU>
SAVE_FIGURES=FIG_POLICY.save_any;
root=fileparts(mfilename('fullpath')); if isempty(root); root=pwd; end
addpath(root); addpath(fullfile(root,'function'));
if exist('cvx_begin','file')~=2; error('CVX is required. Run cvx_setup first.'); end

ITER=cfg.ITER_EXP9;
nmseTarget_dB=cfg.csi_nmse_dB_range_exp9(:).';
epsRange=cfg.csi_error_range_exp9(:).';
assert(numel(nmseTarget_dB)==numel(epsRange),'Exp.9 NMSE/epsilon grids must have equal length.');
E=numel(epsRange);
gamma=10^(cfg.gamma_exp9_dB/10); P_total=cfg.P_total_exp9;
labels={'Fixed RET, No RIS','Fixed RET, RIS-RadarSNR','MA-RCG-MF'};
S=numel(labels);
ObjActual=nan(ITER,E,S); Pcomm=nan(ITER,E,S); MinSINR_dB=nan(ITER,E,S);
DesignFeasible=false(ITER,E,S); QoSSatisfied=false(ITER,E,S); NMSE=nan(ITER,E,5);

fprintf('\nExperiment 9  target-NMSE imperfect CSI (%s mode, MC=%d)\n',upper(RUN_MODE),ITER);
fprintf('Operating point: gamma=%.1f dB, P_T=%.1f dBm\n',cfg.gamma_exp9_dB,cfg.P_total_exp9_dBm);
fprintf('Target NMSE grid (dB): %s\n',strjoin(csi_tick_labels(nmseTarget_dB),', '));
fprintf('Equivalent epsilon:   %s\n',mat2str(epsRange,4));
use_parallel=setup_parallel(cfg,ITER,root);
if use_parallel
    parfor mc=1:ITER
        [o,p,m,fd,q,nm]=run_one(mc,cfg,epsRange,gamma,P_total);
        ObjActual(mc,:,:)=o; Pcomm(mc,:,:)=p; MinSINR_dB(mc,:,:)=m;
        DesignFeasible(mc,:,:)=fd; QoSSatisfied(mc,:,:)=q; NMSE(mc,:,:)=nm;
    end
else
    for mc=1:ITER
        [o,p,m,fd,q,nm]=run_one(mc,cfg,epsRange,gamma,P_total);
        ObjActual(mc,:,:)=o; Pcomm(mc,:,:)=p; MinSINR_dB(mc,:,:)=m;
        DesignFeasible(mc,:,:)=fd; QoSSatisfied(mc,:,:)=q; NMSE(mc,:,:)=nm;
        if SHOW_PROGRESS; fprintf('Exp9 %d/%d\n',mc,ITER); end
    end
end

% -------------------------------------------------------------------------
% Measured NMSE verification. Average in the linear domain across the five
% modeled links and all Monte-Carlo realizations, then convert to dB.
% -------------------------------------------------------------------------
MeasuredNMSE_linear=nan(1,E);
MeasuredNMSE_dB=nan(1,E);
for ee=1:E
    tmp=reshape(NMSE(:,ee,:),[],1);
    tmp=tmp(isfinite(tmp));
    if isempty(tmp)
        continue;
    end
    MeasuredNMSE_linear(ee)=mean(tmp);
    if MeasuredNMSE_linear(ee)==0
        MeasuredNMSE_dB(ee)=-Inf;
    else
        MeasuredNMSE_dB(ee)=10*log10(MeasuredNMSE_linear(ee));
    end
end

% -------------------------------------------------------------------------
% Fixed perfect-CSI cohort avoids survivor bias as the error level grows.
% A sample belongs to the cohort only if ALL compared designs satisfy actual
% QoS at the perfect-CSI point.  The same samples are used at every NMSE.
% -------------------------------------------------------------------------
FixedCohort=all(QoSSatisfied(:,1,:),3);
FixedCohortGrid=repmat(FixedCohort,1,E);
if ~any(FixedCohort)
    warning('Exp9:FixedCohortEmpty','The perfect-CSI common cohort is empty.');
end
if isfield(cfg,'exp9_min_fixed_cohort_warning') && nnz(FixedCohort)<cfg.exp9_min_fixed_cohort_warning
    warning('Exp9:SmallFixedCohort', ...
        'Only %d/%d samples are in the fixed perfect-CSI cohort. Interpret CIs cautiously.', ...
        nnz(FixedCohort),ITER);
end

% Outage-penalized sensing utility: zero contribution whenever actual QoS is
% violated (or the estimated-CSI design is invalid), but only for members of
% the fixed perfect-CSI cohort. Samples outside that cohort remain NaN.
ObjOutagePenalized=nan(size(ObjActual));
for ee=1:E
    for ss=1:S
        v=zeros(ITER,1);
        good=FixedCohort & QoSSatisfied(:,ee,ss) & isfinite(ObjActual(:,ee,ss));
        v(good)=ObjActual(good,ee,ss);
        v(~FixedCohort)=NaN;
        ObjOutagePenalized(:,ee,ss)=v;
    end
end
[muObj,loObj,hiObj,nObj]=bootstrap_mean_ci_to_db(ObjOutagePenalized, ...
    FixedCohortGrid,1,cfg.exp9_bootstrap_reps,cfg.seed_base+990000);

% Conditional sensing statistics are diagnostic only.
CommonQoS=all(QoSSatisfied,3);
[muObjConditional,loObjConditional,hiObjConditional,nObjConditional]= ...
    masked_mean_ci_to_db(ObjActual,CommonQoS,1);

% Reliability is also reported on the SAME fixed perfect-CSI cohort.  The
% all-realization rates are saved separately for diagnostics.
QoSRateAll=100*squeeze(mean(QoSSatisfied,1));
DesignRateAll=100*squeeze(mean(DesignFeasible,1));
QoSRate=nan(E,S); DesignRate=nan(E,S);
for ee=1:E
    for ss=1:S
        if any(FixedCohort)
            QoSRate(ee,ss)=100*mean(QoSSatisfied(FixedCohort,ee,ss));
            DesignRate(ee,ss)=100*mean(DesignFeasible(FixedCohort,ee,ss));
        end
    end
end
[minMu,minLo,minHi]=masked_mean_ci_linear(MinSINR_dB,FixedCohortGrid);
MinSINRMargin_dB=MinSINR_dB-cfg.gamma_exp9_dB;
[marginMu,marginLo,marginHi]=masked_mean_ci_linear(MinSINRMargin_dB,FixedCohortGrid);

% Paired sensing loss is conditional on both perfect and impaired designs
% satisfying QoS for the same sample. It is supplementary, not the main curve.
Loss=nan(ITER,E,S);
for ee=1:E
    for ss=1:S
        ok=FixedCohort & QoSSatisfied(:,ee,ss) & QoSSatisfied(:,1,ss) & ...
            ObjActual(:,ee,ss)>0 & ObjActual(:,1,ss)>0;
        Loss(ok,ee,ss)=10*log10(ObjActual(ok,1,ss)./ObjActual(ok,ee,ss));
    end
end
[lossMu,lossLo,lossHi]=masked_mean_ci_linear(Loss,isfinite(Loss));

% -------------------------------------------------------------------------
% Optional per-experiment figure. The final manuscript uses the composite F4.
% -------------------------------------------------------------------------
if SAVE_MAIN_FIGURE
    run_stamp=datestr(now,'yyyymmdd_HHMMSS');
    figdir=fullfile(root,sprintf('Exp9_Figures_%s_%s',lower(RUN_MODE),run_stamp));
    if SAVE_FIGURES && ~exist(figdir,'dir'); mkdir(figdir); end
    xCSI=1:E; tickLabels=csi_tick_labels(nmseTarget_dB);
    fig=figure('Color','w','Position',[50 50 1500 900]);
    tl=tiledlayout(2,2,'TileSpacing','compact');
    nexttile; hold on;
    for ss=1:S
        errorbar(xCSI,muObj(:,ss),muObj(:,ss)-loObj(:,ss),hiObj(:,ss)-muObj(:,ss), ...
            '-o','LineWidth',1.6);
    end
    configure_nmse_axis(gca,xCSI,tickLabels);
    grid on; ylabel('Outage-penalized weighted radar SNR (dB)');
    legend(labels,'Location','best'); title('(a) Fixed-cohort sensing robustness');
    nexttile; hold on;
    for ss=1:S
        plot(xCSI,QoSRate(:,ss),'-s','LineWidth',1.6,'MarkerSize',7);
    end
    configure_nmse_axis(gca,xCSI,tickLabels);
    grid on; ylim([0 105]); ylabel('Actual QoS satisfaction rate (%)');
    title('(b) Fixed-cohort communication reliability');
    nexttile; hold on;
    for ss=1:S
        errorbar(xCSI,lossMu(:,ss),lossMu(:,ss)-lossLo(:,ss),lossHi(:,ss)-lossMu(:,ss), ...
            '-d','LineWidth',1.5);
    end
    yline(0,'--k'); configure_nmse_axis(gca,xCSI,tickLabels);
    grid on; ylabel('Radar-SNR loss from perfect CSI (dB)');
    title('(c) Paired sensing degradation');
    nexttile; hold on;
    for ss=1:S
        errorbar(xCSI,minMu(:,ss),minMu(:,ss)-minLo(:,ss),minHi(:,ss)-minMu(:,ss), ...
            '-^','LineWidth',1.5);
    end
    yline(cfg.gamma_exp9_dB,'--k','Required QoS');
    configure_nmse_axis(gca,xCSI,tickLabels);
    grid on; ylabel('Minimum actual user SINR (dB)');
    title('(d) Worst-user SINR on fixed cohort');
    sgtitle(tl,sprintf('Imperfect-CSI Sensitivity, QoS=%.1f dB, P_T=%.1f dBm', ...
        cfg.gamma_exp9_dB,cfg.P_total_exp9_dBm));
    save_figure(fig,figdir,'01_CSI_NMSE_Sensitivity');
end

% -------------------------------------------------------------------------
% Save raw data and a publication-ready table.
% -------------------------------------------------------------------------
results_dir=fullfile(root,'results'); if ~exist(results_dir,'dir'); mkdir(results_dir); end
T=table(); rowCount=0;
for ee=1:E
    for ss=1:S
        row=table(nmseTarget_dB(ee),epsRange(ee),MeasuredNMSE_dB(ee),string(labels{ss}), ...
            muObj(ee,ss),loObj(ee,ss),hiObj(ee,ss),muObjConditional(ee,ss), ...
            QoSRate(ee,ss),QoSRateAll(ee,ss),DesignRate(ee,ss),DesignRateAll(ee,ss), ...
            minMu(ee,ss),marginMu(ee,ss),lossMu(ee,ss),nObj(ee,ss),nObjConditional(ee,ss), ...
            'VariableNames',{'Target_NMSE_dB','ErrorCoefficient_epsilon','MeasuredNMSE_Mean_dB', ...
            'Scheme','OutagePenalizedRadarSNR_Mean_dB','CI95_Low_dB','CI95_High_dB', ...
            'ConditionalRadarSNR_Mean_dB','FixedCohort_ActualQoS_percent', ...
            'AllSamples_ActualQoS_percent','FixedCohort_DesignFeasible_percent', ...
            'AllSamples_DesignFeasible_percent','MinSINR_Mean_dB','MinSINR_Margin_Mean_dB', ...
            'ConditionalRadarLoss_Mean_dB','FixedCohort_samples','ConditionalQoS_samples'});
        rowCount=rowCount+1;
        if rowCount==1; T=row; else; T=[T;row]; end %#ok<AGROW>
    end
end
matfile=fullfile(results_dir,sprintf('exp9_csi_uncertainty_%s.mat',lower(RUN_MODE)));
save(matfile,'cfg','nmseTarget_dB','epsRange','labels','ObjActual','ObjOutagePenalized', ...
    'Pcomm','MinSINR_dB','DesignFeasible','QoSSatisfied','NMSE', ...
    'MeasuredNMSE_linear','MeasuredNMSE_dB','FixedCohort','FixedCohortGrid','CommonQoS', ...
    'muObj','loObj','hiObj','nObj','muObjConditional','loObjConditional', ...
    'hiObjConditional','nObjConditional','QoSRate','QoSRateAll','DesignRate','DesignRateAll', ...
    'Loss','lossMu','lossLo','lossHi','minMu','minLo','minHi', ...
    'MinSINRMargin_dB','marginMu','marginLo','marginHi','T');
writetable(T,fullfile(results_dir,sprintf('exp9_csi_uncertainty_%s.csv',lower(RUN_MODE))));

fprintf('Exp.9 fixed perfect-CSI cohort: %d/%d samples.\n',nnz(FixedCohort),ITER);
fprintf('Measured mean NMSE (dB): %s\n',strjoin(csi_tick_labels(MeasuredNMSE_dB),', '));
fprintf('Fixed-cohort QoS at strongest error (%%): %s\n',mat2str(QoSRate(end,:),4));
fprintf('Experiment 9 complete in %.1f s. Results: %s\n',toc(t_total),matfile);

function [Obj,Pc,MinSINR,DesignFea,QoS,NMSE]=run_one(mc,cfg,epsRange,gamma,P_total)
E=numel(epsRange); S=3; Obj=nan(1,E,S); Pc=Obj; MinSINR=Obj;
DesignFea=false(1,E,S); QoS=false(1,E,S); NMSE=nan(1,E,5);
seed=cfg.seed_base+cfg.seed_stride*9+mc;
scenTrue=generate_isac_scenario(cfg,cfg.N,cfg.K,seed,[],[],struct());
for ee=1:E
    [scenEst,rep]=perturb_scenario_csi(scenTrue,epsRange(ee),seed+900000);
    NMSE(1,ee,:)=[rep.nmse_Hd_user rep.nmse_Hr_user rep.nmse_Hd_tar rep.nmse_Hr_tar rep.nmse_G];

    % 1) Fixed RET, no RIS.
    r1=evaluate_fixed_ret_ris_config(cfg,scenEst,cfg.theta_ref,[],gamma,P_total,false);
    [Obj(1,ee,1),Pc(1,ee,1),MinSINR(1,ee,1),DesignFea(1,ee,1),QoS(1,ee,1)]= ...
        test_design(cfg,scenTrue,r1,[],false,gamma);

    % 2) Fixed RET, radar-oriented RIS.
    cfix=cfg; cfix.exp6_fresh_candidate_mode='full'; cfix.exp6_fresh_allow_random=true;
    rng(seed+910000,'twister');
    [phi2,r2]=optimize_fixed_ret_ris_radar(cfix,scenEst,cfg.theta_ref,ones(cfg.N,1),gamma,P_total);
    [Obj(1,ee,2),Pc(1,ee,2),MinSINR(1,ee,2),DesignFea(1,ee,2),QoS(1,ee,2)]= ...
        test_design(cfg,scenTrue,r2,phi2,true,gamma);

    % 3) Joint RET-RIS radar-oriented design.
    rng(seed+920000,'twister'); out=evaluate_ret_ris_schemes(cfg,scenEst,gamma,P_total);
    if out.Feasible(5)
        r3=evaluate_fixed_ret_ris_config(cfg,scenEst,out.Theta(5),out.Phi{5},gamma,P_total,true);
    else
        r3=struct('feasible',false);
    end
    [Obj(1,ee,3),Pc(1,ee,3),MinSINR(1,ee,3),DesignFea(1,ee,3),QoS(1,ee,3)]= ...
        test_design(cfg,scenTrue,r3,get_phi(out,5),true,gamma);
end
end

function [obj,pc,minsinr,designfea,qos]=test_design(cfg,scenTrue,r,phi,use_ris,gamma)
obj=NaN; pc=NaN; minsinr=NaN; designfea=false; qos=false;
if ~isstruct(r)||~isfield(r,'feasible')||~r.feasible; return; end
designfea=true;
a=evaluate_external_beamformers(cfg,scenTrue,r.theta_deg,phi,use_ris,r.Wc,r.Wr,gamma);
if ~a.valid; return; end
obj=a.obj_snr; pc=a.p_comm; minsinr=a.min_sinr_dB; qos=a.qos_satisfied;
end

function phi=get_phi(out,s)
phi=[]; if isstruct(out)&&isfield(out,'Phi')&&numel(out.Phi)>=s; phi=out.Phi{s}; end
end

function labels=csi_tick_labels(v)
labels=cell(1,numel(v));
for ii=1:numel(v)
    if isinf(v(ii)) && v(ii)<0
        labels{ii}='Perfect';
    elseif isfinite(v(ii))
        labels{ii}=sprintf('%g',v(ii));
    else
        labels{ii}='NaN';
    end
end
end

function configure_nmse_axis(ax,x,tickLabels)
set(ax,'XTick',x,'XTickLabel',tickLabels);
xlim(ax,[0.7 numel(x)+0.3]);
xlabel(ax,'Target CSI NMSE (dB)');
end

function use_parallel=setup_parallel(cfg,ITER,root)
use_parallel=false;
try
    pool=gcp('nocreate'); wanted=min(cfg.num_workers,ITER);
    if isempty(pool)
        pool=parpool('Processes',wanted);
    elseif pool.NumWorkers~=wanted
        delete(pool); pool=parpool('Processes',wanted);
    end
    addAttachedFiles(pool,{fullfile(root,'unified_ret_ris_isac_config.m')});
    use_parallel=true;
catch ME
    warning('Exp9:ParallelFallback','Parallel pool unavailable (%s). Running serially.',ME.message);
end
end

function save_figure(fig,dirpath,name)
if ~exist(dirpath,'dir'); mkdir(dirpath); end
try; exportgraphics(fig,fullfile(dirpath,[name '.pdf']),'ContentType','vector'); catch; end
try; exportgraphics(fig,fullfile(dirpath,[name '.png']),'Resolution',300); catch; end
try; savefig(fig,fullfile(dirpath,[name '.fig'])); catch; end
end
