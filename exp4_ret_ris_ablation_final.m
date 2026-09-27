% =========================================================================
% Experiment 4 FINAL: RET-RIS Dual-Domain Ablation Study
%
% Purpose:
%   Quantify the separate and joint contributions of
%     1) communication-power release, and
%     2) sensing-channel shaping.
%
% A representative operating point is used:
%   gamma = cfg.gamma_exp4_dB, P_T = cfg.P_total_exp4_dBm.
%
% Five schemes:
%   1) Fixed RET, No RIS
%   2) Fixed RET, RIS-MinPc
%   3) RET-NoRIS-MinPc
%   4) RET-RIS-MinPc
%   5) RET-RIS-RadarSNR
%
% Six paper figures are generated:
%   01 communication power
%   02 residual radar-power share
%   03 unit-power sensing gain
%   04 final weighted radar SNR
%   05 communication/radar waveform sensing contribution
%   06 dual-mechanism map
% =========================================================================

clear; close all; clc;

%% ======================== Controls =======================================
env_mode=getenv('RET_RIS_RUN_MODE');
if isempty(env_mode); RUN_MODE='paper'; else; RUN_MODE=lower(env_mode); end
SHOW_PROGRESS = true;
t_total = tic;

cfg = unified_ret_ris_isac_config(RUN_MODE);
FIG_POLICY=get_figure_policy(RUN_MODE);
SAVE_MAIN_FIGURE=FIG_POLICY.save_main;
SAVE_SUPP_FIGURES=FIG_POLICY.save_supplementary;
SAVE_FIGURES=FIG_POLICY.save_any;
ITER = cfg.ITER_EXP4;
gamma_dB = cfg.gamma_exp4_dB;
gamma = 10^(gamma_dB/10);
P_total_dBm = cfg.P_total_exp4_dBm;
P_total = cfg.P_total_exp4;
num_workers = min(cfg.num_workers,ITER);

%% ======================== Paths / Checks =================================
root=fileparts(mfilename('fullpath'));
if isempty(root); root=pwd; end
addpath(root);
func_path=fullfile(root,'function');
if exist(func_path,'dir'); addpath(func_path); end

required={'generate_isac_scenario','evaluate_ret_ris_schemes', ...
    'generate_channel','solve_comm_min_power','build_radar_matrix', ...
    'build_effective_user_channel','allocate_radar_power', ...
    'optimize_ret_only_radar','optimize_fixed_ret_ris_radar','opt_phi_min_pc','opt_phi_radar_gain','opt_a_min_pc'};
for ii=1:numel(required)
    if exist(required{ii},'file')~=2
        error('Required function not found: %s.m',required{ii});
    end
end
if exist('cvx_begin','file')~=2
    error('CVX is required. Install CVX and run cvx_setup first.');
end

stamp=datestr(now,'yyyymmdd_HHMMSS');
figure_dir=fullfile(root,sprintf('Exp4_Figures_%s_%s',lower(RUN_MODE),stamp));
if SAVE_FIGURES && ~exist(figure_dir,'dir'); mkdir(figure_dir); end

%% ======================== Parallel Pool ==================================
use_parallel=true;
if use_parallel
    try
        poolobj=gcp('nocreate');
        if ~isempty(poolobj) && poolobj.NumWorkers~=num_workers
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

fprintf('\n===============================================================\n');
fprintf(' EXPERIMENT 4: RET-RIS DUAL-DOMAIN ABLATION\n');
fprintf(' Mode=%s, MC=%d, gamma=%.1f dB, P_T=%.1f dBm\n', ...
    upper(RUN_MODE),ITER,gamma_dB,P_total_dBm);
fprintf(' M=%d, N=%d, K=%d, T=%d\n',cfg.M,cfg.N,cfg.K,cfg.T);
if SAVE_FIGURES; fprintf(' Figures: %s\n',figure_dir); end
fprintf('===============================================================\n\n');

%% ======================== Buffers ========================================
S=5;
Pcomm=nan(ITER,S);
Pradar=nan(ITER,S);
Obj=nan(ITER,S);
G_eff=nan(ITER,S);
Jcomm=nan(ITER,S);
Jradar=nan(ITER,S);
Theta=nan(ITER,S);
Feasible=false(ITER,S);
Runtime=nan(ITER,S);

% Same-objective hardware ablation: all four schemes maximize radar SNR.
Ssame=4;
ObjSame=nan(ITER,Ssame); PcommSame=nan(ITER,Ssame); PradarSame=nan(ITER,Ssame);
GeffSame=nan(ITER,Ssame); FeasibleSame=false(ITER,Ssame);
ThetaSame=nan(ITER,Ssame); RuntimeSame=nan(ITER,Ssame);

%% ======================== Progress =======================================
progressQueue=[];
if SHOW_PROGRESS
    try
        progressQueue=parallel.pool.DataQueue;
        afterEach(progressQueue,@exp4_progress_bar);
        exp4_progress_bar(struct('action','reset','total',ITER));
    catch ME
        warning('Progress bar initialization failed: %s',ME.message);
        progressQueue=[];
    end
end

%% ======================== Monte Carlo ====================================
t_mc=tic;
if use_parallel
    parfor iter=1:ITER
        seed=cfg.seed_base+iter*cfg.seed_stride;
        scen=generate_isac_scenario(cfg,cfg.N,cfg.K,seed,[],[],struct());
        r=evaluate_ret_ris_schemes(cfg,scen,gamma,P_total);
        Pcomm(iter,:)=r.Pcomm;
        Pradar(iter,:)=r.Pradar;
        Obj(iter,:)=r.Obj;
        G_eff(iter,:)=r.G_eff;
        Jcomm(iter,:)=r.Jcomm;
        Jradar(iter,:)=r.Jradar;
        Theta(iter,:)=r.Theta;
        Feasible(iter,:)=r.Feasible;
        Runtime(iter,:)=r.RuntimeCumulative;
        [ObjSame(iter,:),PcommSame(iter,:),PradarSame(iter,:), ...
            GeffSame(iter,:),FeasibleSame(iter,:),ThetaSame(iter,:), ...
            RuntimeSame(iter,:)]=run_same_objective_ablation(cfg,scen,gamma,P_total,r,seed);
        if ~isempty(progressQueue); send(progressQueue,1); end
    end
else
    for iter=1:ITER
        seed=cfg.seed_base+iter*cfg.seed_stride;
        scen=generate_isac_scenario(cfg,cfg.N,cfg.K,seed,[],[],struct());
        r=evaluate_ret_ris_schemes(cfg,scen,gamma,P_total);
        Pcomm(iter,:)=r.Pcomm;
        Pradar(iter,:)=r.Pradar;
        Obj(iter,:)=r.Obj;
        G_eff(iter,:)=r.G_eff;
        Jcomm(iter,:)=r.Jcomm;
        Jradar(iter,:)=r.Jradar;
        Theta(iter,:)=r.Theta;
        Feasible(iter,:)=r.Feasible;
        Runtime(iter,:)=r.RuntimeCumulative;
        [ObjSame(iter,:),PcommSame(iter,:),PradarSame(iter,:), ...
            GeffSame(iter,:),FeasibleSame(iter,:),ThetaSame(iter,:), ...
            RuntimeSame(iter,:)]=run_same_objective_ablation(cfg,scen,gamma,P_total,r,seed);
        if ~isempty(progressQueue); send(progressQueue,1); end
    end
end
runtime=toc(t_mc);

%% ======================== Statistics =====================================
labels={'Fixed RET, No RIS','Fixed RET, RIS-MinPc', ...
        'RET-NoRIS-MinPc','RET-RIS-MinPc','RET-RIS-RadarSNR'};
short_labels={'Fixed/No RIS','RIS only','RET only','RET+RIS MinPc','RET+RIS RadarSNR'};
% Keep non-plot metadata outside figure-policy branches. PAPER master mode
% suppresses per-experiment figures, but these labels are still required by
% the same-objective ablation table and MAT/CSV export below.
same_labels={'Fixed RET/No RIS','RET-only RadarSNR', ...
    'Fixed RET/RIS-RadarSNR','Joint RET-RIS RadarSNR'};

MechanismCommon=all(Feasible,2);
MechanismCommonRate=100*mean(MechanismCommon);
[Pcomm_dBm,Pcomm_lo,Pcomm_hi,~]=masked_mean_ci_to_db( ...
    reshape(Pcomm,[ITER 1 S]),MechanismCommon,1e-3);
Pcomm_dBm=Pcomm_dBm.'; Pcomm_lo=Pcomm_lo.'; Pcomm_hi=Pcomm_hi.';
[Pradar_pct,Pradar_lo,Pradar_hi,~]=masked_mean_ci_linear( ...
    reshape(100*Pradar/P_total,[ITER 1 S]),MechanismCommon);
Pradar_pct=Pradar_pct.'; Pradar_lo=Pradar_lo.'; Pradar_hi=Pradar_hi.';
[Geff_dB,Geff_lo,Geff_hi,~]=masked_mean_ci_to_db( ...
    reshape(G_eff,[ITER 1 S]),MechanismCommon,1);
Geff_dB=Geff_dB.'; Geff_lo=Geff_lo.'; Geff_hi=Geff_hi.';
[Obj_dB,Obj_lo,Obj_hi,~]=masked_mean_ci_to_db( ...
    reshape(Obj,[ITER 1 S]),MechanismCommon,1);
Obj_dB=Obj_dB.'; Obj_lo=Obj_lo.'; Obj_hi=Obj_hi.';

Fea_rate=100*mean(Feasible,1);
Runtime_avg=nanmean_vector(Runtime);

SameCommon=all(FeasibleSame,2);
SameCommonRate=100*mean(SameCommon);
[ObjSame_dB,ObjSame_lo,ObjSame_hi,ObjSame_n]= ...
    masked_mean_ci_to_db(reshape(ObjSame,[ITER 1 Ssame]),SameCommon,1);
ObjSame_dB=ObjSame_dB.'; ObjSame_lo=ObjSame_lo.'; ObjSame_hi=ObjSame_hi.'; ObjSame_n=ObjSame_n.';
[PcommSame_dBm,PcommSame_lo,PcommSame_hi,~]= ...
    masked_mean_ci_to_db(reshape(PcommSame,[ITER 1 Ssame]),SameCommon,1e-3);
PcommSame_dBm=PcommSame_dBm.'; PcommSame_lo=PcommSame_lo.'; PcommSame_hi=PcommSame_hi.';
[GeffSame_dB,GeffSame_lo,GeffSame_hi,~]= ...
    masked_mean_ci_to_db(reshape(GeffSame,[ITER 1 Ssame]),SameCommon,1);
GeffSame_dB=GeffSame_dB.'; GeffSame_lo=GeffSame_lo.'; GeffSame_hi=GeffSame_hi.';
SameFeaRate=100*mean(FeasibleSame,1);
SameRuntimeMed=nanmedian_vector(RuntimeSame);

contrib_comm=nan(1,S);
contrib_radar=nan(1,S);
for s=1:S
    total=Jcomm(:,s)+Jradar(:,s);
    valid=MechanismCommon & isfinite(total) & total>0;
    if any(valid)
        contrib_comm(s)=100*mean(Jcomm(valid,s)./total(valid));
        contrib_radar(s)=100*mean(Jradar(valid,s)./total(valid));
    end
end

if SAVE_MAIN_FIGURE
%% ======================== Main Paper Composite ===========================
% The main paper uses one 2x2 mechanism figure. Residual-power and waveform
% contribution plots are retained as supplementary diagnostics.
fig0=figure('Position',[30 30 1450 900],'Color','w');
tl0=tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

nexttile;
bar_with_ci(Pcomm_dBm,Pcomm_lo,Pcomm_hi,short_labels); grid on;
ylabel('Communication power P_c^{min} (dBm)'); title('(a) Communication-power release');

nexttile;
bar_with_ci(Geff_dB,Geff_lo,Geff_hi,short_labels); grid on;
ylabel('Unit-power sensing gain G_{eff} (dB/W)'); title('(b) Sensing-channel shaping');

nexttile;
bar_with_ci(Obj_dB,Obj_lo,Obj_hi,short_labels); grid on;
ylabel('Weighted radar SNR (dB)'); title('(c) Final sensing performance');

nexttile; hold on;
mk={'o','^','d','s','p'};
for s=1:S
    scatter(Pradar_pct(s),Geff_dB(s),120,mk{s},'LineWidth',1.7);
    text(Pradar_pct(s)+0.35,Geff_dB(s),short_labels{s}, ...
        'FontSize',9,'VerticalAlignment','middle');
end
grid on; box on; xlabel('Residual sensing-power share 100P_r/P_T (%)');
ylabel('Unit-power sensing gain G_{eff} (dB/W)');
title('(d) Dual-mechanism map');

sgtitle(tl0,sprintf('Mechanism Analysis of RET-RIS Cooperation (QoS threshold = %.1f dB)',gamma_dB));
if SAVE_MAIN_FIGURE; save_figure_auto(fig0,figure_dir,'00_Main_Dual_Mechanism_Composite'); end

end

if SAVE_SUPP_FIGURES
%% ======================== Figure 1: Communication Power ==================
fig1=figure('Position',[60 60 1050 620],'Color','w');
bar_with_ci(Pcomm_dBm,Pcomm_lo,Pcomm_hi,short_labels);
ylabel('Minimum communication power P_c^{min} (dBm)');
title(sprintf('Ablation: Communication-Power Release (QoS = %.1f dB)',gamma_dB));
grid on;
if SAVE_SUPP_FIGURES; save_figure_auto(fig1,figure_dir,'S1_Communication_Power_Ablation'); end

%% ======================== Figure 2: Residual Radar Power =================
fig2=figure('Position',[80 80 1050 620],'Color','w');
bar_with_ci(Pradar_pct,Pradar_lo,Pradar_hi,short_labels);
ylabel('Residual radar-power share 100P_r/P_T (%)');
title('Ablation: Power Released to Sensing');
grid on;
if SAVE_SUPP_FIGURES; save_figure_auto(fig2,figure_dir,'S2_Residual_Radar_Power_Ablation'); end

%% ======================== Figure 3: Unit-Power Sensing Gain ===============
fig3=figure('Position',[100 100 1050 620],'Color','w');
bar_with_ci(Geff_dB,Geff_lo,Geff_hi,short_labels);
ylabel('Unit-power sensing gain G_{eff} (dB/W)');
title('Ablation: Sensing-Channel Shaping Gain');
grid on;
if SAVE_SUPP_FIGURES; save_figure_auto(fig3,figure_dir,'S3_Unit_Power_Sensing_Gain'); end

%% ======================== Figure 4: Final Radar SNR =======================
fig4=figure('Position',[120 120 1050 620],'Color','w');
bar_with_ci(Obj_dB,Obj_lo,Obj_hi,short_labels);
ylabel('Weighted radar SNR (dB)');
title('Ablation: Final Sensing Performance');
grid on;
if SAVE_SUPP_FIGURES; save_figure_auto(fig4,figure_dir,'S4_Final_Radar_SNR_Ablation'); end

%% ======================== Figure 5: ISAC Waveform Contribution ============
fig5=figure('Position',[140 140 1050 620],'Color','w');
bar(1:S,[contrib_comm(:),contrib_radar(:)],'stacked');
set(gca,'XTick',1:S,'XTickLabel',short_labels,'XTickLabelRotation',18);
ylabel('Contribution to weighted radar SNR (%)');
legend({'Communication waveforms','Dedicated radar waveform'}, ...
    'Location','northwest');
ylim([0 100]); grid on;
title('Communication and Dedicated-Sensing Contributions to Radar SNR');
if SAVE_SUPP_FIGURES; save_figure_auto(fig5,figure_dir,'S5_Waveform_Sensing_Contribution'); end

%% ======================== Figure 6: Dual-Mechanism Map ====================
fig6=figure('Position',[160 160 1050 650],'Color','w');
hold on;
mk={'o','^','d','s','p'};
for s=1:S
    scatter(Pradar_pct(s),Geff_dB(s),140,mk{s},'LineWidth',1.8);
    text(Pradar_pct(s)+0.5,Geff_dB(s),short_labels{s}, ...
        'FontSize',10,'VerticalAlignment','middle');
end
xlabel('Residual radar-power share 100P_r/P_T (%)');
ylabel('Unit-power sensing gain G_{eff} (dB/W)');
title('Dual Mechanisms: Power Release versus Sensing-Channel Shaping');
grid on; box on;
if SAVE_SUPP_FIGURES; save_figure_auto(fig6,figure_dir,'S6_Dual_Mechanism_Map'); end

%% ======================== Figure 7: Same-Objective Hardware Ablation =====
fig7=figure('Position',[180 180 1380 620],'Color','w');
tl7=tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
nexttile; bar_with_ci(ObjSame_dB,ObjSame_lo,ObjSame_hi,same_labels); grid on;
ylabel('Weighted radar SNR (dB)');
title(sprintf('Same objective, common n=%d',min(ObjSame_n)));
nexttile; bar_with_ci(PcommSame_dBm,PcommSame_lo,PcommSame_hi,same_labels); grid on;
ylabel('Communication power (dBm)'); title('Communication cost');
nexttile; bar_with_ci(GeffSame_dB,GeffSame_lo,GeffSame_hi,same_labels); grid on;
ylabel('Unit-power sensing gain (dB/W)'); title('Sensing-channel shaping');
sgtitle(tl7,sprintf('Hardware Ablation under the Same Radar-SNR Objective (common rate %.1f%%)', ...
    SameCommonRate));
if SAVE_SUPP_FIGURES; save_figure_auto(fig7,figure_dir,'S7_Same_Objective_Hardware_Ablation'); end


end

%% ======================== Console Tables / Save ===========================
fprintf('\nMechanism ablation summary:\n');
Tsummary=table(labels(:),Pcomm_dBm(:),Pradar_pct(:),Geff_dB(:), ...
    Obj_dB(:),Fea_rate(:),Runtime_avg(:), ...
    'VariableNames',{'Scheme','Pcomm_dBm','Pradar_percent', ...
    'UnitPowerGain_dB','RadarSNR_dB','Feasibility_percent','Runtime_s'});
disp(Tsummary);
Tsame=table(same_labels(:),ObjSame_dB(:),ObjSame_lo(:),ObjSame_hi(:), ...
    PcommSame_dBm(:),GeffSame_dB(:),SameFeaRate(:),SameRuntimeMed(:), ...
    'VariableNames',{'Scheme','RadarSNR_Mean_dB','RadarSNR_CI95_Low_dB', ...
    'RadarSNR_CI95_High_dB','Pcomm_dBm','UnitPowerGain_dB', ...
    'Feasibility_percent','Runtime_median_s'});
fprintf('\nSame-objective hardware ablation (common-feasible rate %.1f%%):\n',SameCommonRate);
disp(Tsame);

results_dir=fullfile(root,'results');
if ~exist(results_dir,'dir'); mkdir(results_dir); end
save(fullfile(results_dir,sprintf('exp4_ablation_revised_%s.mat',lower(RUN_MODE))), ...
    'cfg','gamma_dB','gamma','P_total','P_total_dBm', ...
    'Pcomm','Pradar','Obj','G_eff','Jcomm','Jradar','Theta', ...
    'Feasible','Runtime','MechanismCommon','MechanismCommonRate','Tsummary','ObjSame','PcommSame','PradarSame', ...
    'GeffSame','FeasibleSame','ThetaSame','RuntimeSame','SameCommon', ...
    'SameCommonRate','ObjSame_dB','ObjSame_lo','ObjSame_hi','Tsame');
writetable(Tsummary,fullfile(results_dir,sprintf('exp4_mechanism_summary_%s.csv',lower(RUN_MODE))));
writetable(Tsame,fullfile(results_dir,sprintf('exp4_same_objective_ablation_%s.csv',lower(RUN_MODE))));

total_time=toc(t_total);
fprintf('\n===============================================================\n');
fprintf(' EXPERIMENT 4 COMPLETED\n');
fprintf(' Monte Carlo time: %.1f s (%.2f min)\n',runtime,runtime/60);
fprintf(' Total runtime:     %.1f s (%.2f min)\n',total_time,total_time/60);
fprintf(' Average per MC:    %.1f s\n',runtime/ITER);
if SAVE_FIGURES; fprintf(' Figure folder:     %s\n',figure_dir); end
fprintf('===============================================================\n');

%% ========================================================================
%% Local helpers
%% ========================================================================

function [obj,pc,pr,geff,fea,theta,trun]=run_same_objective_ablation(cfg,scen,gamma,P_total,r,seed)
obj=nan(1,4); pc=nan(1,4); pr=nan(1,4); geff=nan(1,4);
fea=false(1,4); theta=nan(1,4); trun=nan(1,4);
obj(1)=r.Obj(1); pc(1)=r.Pcomm(1); pr(1)=r.Pradar(1); geff(1)=r.G_eff(1);
fea(1)=r.Feasible(1); theta(1)=r.Theta(1); trun(1)=r.RuntimeCumulative(1);
t=tic; [theta2,res2]=optimize_ret_only_radar(cfg,scen,gamma,P_total); trun(2)=toc(t);
[obj(2),pc(2),pr(2),geff(2),fea(2)]=extract_result_metrics(res2); theta(2)=theta2;
cfg_fixed=cfg;
cfg_fixed.exp6_fresh_candidate_mode='full';
cfg_fixed.exp6_fresh_allow_random=true;
cfg_fixed.exp6_fresh_ris_trials=cfg.num_phase_trials;
rng(seed+440000,'twister');
t=tic; [~,res3]=optimize_fixed_ret_ris_radar(cfg_fixed,scen,cfg.theta_ref, ...
    ones(scen.N,1),gamma,P_total); trun(3)=toc(t);
[obj(3),pc(3),pr(3),geff(3),fea(3)]=extract_result_metrics(res3); theta(3)=cfg.theta_ref;
obj(4)=r.Obj(5); pc(4)=r.Pcomm(5); pr(4)=r.Pradar(5); geff(4)=r.G_eff(5);
fea(4)=r.Feasible(5); theta(4)=r.Theta(5); trun(4)=r.RuntimeCumulative(5);
end

function [obj,pc,pr,geff,fea]=extract_result_metrics(res)
obj=NaN; pc=NaN; pr=NaN; geff=NaN; fea=false;
if ~isstruct(res)||~isfield(res,'feasible')||~res.feasible; return; end
fea=true; obj=res.obj_snr; pc=res.p_comm; pr=res.p_radar;
Z=null(res.Hu');
if isempty(Z); geff=0; return; end
B=(Z'*res.C*Z); B=(B+B')/2; geff=max(max(real(eig(B))),0);
end

function m=nanmedian_vector(X)
m=nan(1,size(X,2));
for jj=1:size(X,2)
    v=X(:,jj); v=v(isfinite(v)); if ~isempty(v); m(jj)=median(v); end
end
end

function [mu_db,lo_db,hi_db]=linear_ci_to_db(X,reference)
[mu,lo,hi]=linear_ci(X);
mu_db=nan(size(mu)); lo_db=mu_db; hi_db=mu_db;
v=isfinite(mu)&mu>0;
mu_db(v)=10*log10(mu(v)/reference);
lo_db(v)=10*log10(max(lo(v),eps)/reference);
hi_db(v)=10*log10(max(hi(v),eps)/reference);
end

function [mu,lo,hi]=linear_ci(X)
S=size(X,2);
mu=nan(1,S); lo=mu; hi=mu;
for s=1:S
    v=X(:,s); v=v(isfinite(v));
    if isempty(v); continue; end
    mu(s)=mean(v);
    if numel(v)>1
        half=student_t_critical(0.95,numel(v)-1)*std(v,0)/sqrt(numel(v));
    else
        half=0;
    end
    lo(s)=mu(s)-half;
    hi(s)=mu(s)+half;
end
end

function m=nanmean_vector(X)
m=nan(1,size(X,2));
for s=1:size(X,2)
    v=X(:,s); v=v(isfinite(v));
    if ~isempty(v); m(s)=mean(v); end
end
end

function bar_with_ci(mu,lo,hi,xlabels)
S=numel(mu);
bar(1:S,mu); hold on;
elow=mu-lo; ehigh=hi-mu;
errorbar(1:S,mu,elow,ehigh,'k','LineStyle','none','LineWidth',1.2);
set(gca,'XTick',1:S,'XTickLabel',xlabels,'XTickLabelRotation',18);
xlim([0.4 S+0.6]);
end

function save_figure_auto(fig,folder,name)
if ~exist(folder,'dir'); mkdir(folder); end
set(fig,'Color','w'); drawnow;
png=fullfile(folder,[name '.png']);
figfile=fullfile(folder,[name '.fig']);
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

function exp4_progress_bar(message)
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
        if completed>0
            eta=elapsed/completed*(total-completed);
        else
            eta=NaN;
        end
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
