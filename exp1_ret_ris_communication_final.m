% =========================================================================
% Experiment 1 REVISED: Communication power saving with RET and RIS
%
% Four communication-oriented schemes are evaluated on identical channel
% realizations. All plotted power means use the samples for which all four
% schemes are feasible. Marginal gains are computed realization by
% realization in dB, and 95% confidence intervals are reported.
% =========================================================================
clear; close all; clc;
t_total=tic; SHOW_PROGRESS=true;

env_mode=getenv('RET_RIS_RUN_MODE');
if isempty(env_mode); RUN_MODE='paper'; else; RUN_MODE=lower(env_mode); end
cfg=unified_ret_ris_isac_config(RUN_MODE);
FIG_POLICY=get_figure_policy(RUN_MODE);
SAVE_MAIN_FIGURE=FIG_POLICY.save_main;
SAVE_SUPP_FIGURES=FIG_POLICY.save_supplementary;
SAVE_FIGURES=FIG_POLICY.save_any;
ITER=cfg.ITER; gamma_dB_range=cfg.gamma_dB_range_exp1;
gamma_range=10.^(gamma_dB_range/10); numG=numel(gamma_range);
S=4; labels={'Fixed RET, No RIS','Fixed RET, RIS-MinPc', ...
    'RET-MinPc, No RIS','RET-RIS-MinPc'};
markers={'-o','-^','-d','-s'};

root=fileparts(mfilename('fullpath')); if isempty(root); root=pwd; end
addpath(root); addpath(fullfile(root,'function'));
if exist('cvx_begin','file')~=2
    error('CVX is required. Install CVX and run cvx_setup first.');
end
stamp=datestr(now,'yyyymmdd_HHMMSS');
figure_dir=fullfile(root,sprintf('Exp1_Figures_%s_%s',lower(RUN_MODE),stamp));
if SAVE_FIGURES && ~exist(figure_dir,'dir'); mkdir(figure_dir); end

fprintf('\n===============================================================\n');
fprintf(' EXPERIMENT 1: PAIRED COMMUNICATION-POWER SAVING\n');
fprintf(' Mode=%s, MC=%d, M=%d, N=%d, K=%d\n',upper(RUN_MODE), ...
    ITER,cfg.M,cfg.N,cfg.K);
fprintf(' Fixed RET baseline = %.1f deg; hardware grid = %.1f deg\n', ...
    cfg.theta_ref,cfg.ret_angle_resolution_deg);
fprintf('===============================================================\n\n');

Pc=nan(ITER,numG,S); Feasible=false(ITER,numG,S); Theta=nan(ITER,numG,S);
progressQueue=[];
if SHOW_PROGRESS
    try
        progressQueue=parallel.pool.DataQueue;
        afterEach(progressQueue,@(~)exp1_progress(ITER*numG));
    catch; progressQueue=[]; end
end

use_parallel=true; num_workers=min(cfg.num_workers,ITER);
if use_parallel
    try
        poolobj=gcp('nocreate');
        if ~isempty(poolobj) && poolobj.NumWorkers~=num_workers
            delete(poolobj); poolobj=[];
        end
        if isempty(poolobj); parpool(num_workers); end
        pctRunOnAll(['addpath(''' root ''');']);
        pctRunOnAll(['addpath(''' fullfile(root,'function') ''');']);
    catch ME
        warning('Parallel setup failed: %s. Running serially.',ME.message);
        use_parallel=false;
    end
end

t_mc=tic;
if use_parallel
    parfor iter=1:ITER
        seed=cfg.seed_base+iter*cfg.seed_stride;
        scen=generate_isac_scenario(cfg,cfg.N,cfg.K,seed,[],[],struct());
        solver_seed=seed+cfg.solver_seed_offset_sweep;
        sw=evaluate_ret_ris_sweep(cfg,scen,gamma_range,cfg.P_total_exp1, ...
            solver_seed,struct('use_continuation',cfg.exp23_use_continuation, ...
            'seed_policy','operating_point','evaluator_opts',struct('max_scheme',4)));
        pc_i=sw.Pcomm(:,1:S); fea_i=sw.Feasible(:,1:S); th_i=sw.Theta(:,1:S);
        for gg=1:numG
            if ~isempty(progressQueue); send(progressQueue,1); end
        end
        Pc(iter,:,:)=reshape(pc_i,[1 numG S]);
        Feasible(iter,:,:)=reshape(fea_i,[1 numG S]);
        Theta(iter,:,:)=reshape(th_i,[1 numG S]);
    end
else
    for iter=1:ITER
        seed=cfg.seed_base+iter*cfg.seed_stride;
        scen=generate_isac_scenario(cfg,cfg.N,cfg.K,seed,[],[],struct());
        solver_seed=seed+cfg.solver_seed_offset_sweep;
        sw=evaluate_ret_ris_sweep(cfg,scen,gamma_range,cfg.P_total_exp1, ...
            solver_seed,struct('use_continuation',cfg.exp23_use_continuation, ...
            'seed_policy','operating_point','evaluator_opts',struct('max_scheme',4)));
        Pc(iter,:,:)=reshape(sw.Pcomm(:,1:S),[1 numG S]);
        Feasible(iter,:,:)=reshape(sw.Feasible(:,1:S),[1 numG S]);
        Theta(iter,:,:)=reshape(sw.Theta(:,1:S),[1 numG S]);
        for gg=1:numG
            if SHOW_PROGRESS; exp1_progress(ITER*numG); end
        end
    end
end
runtime=toc(t_mc);

%% Paired/common-feasible statistics
CommonFeasible=all(Feasible,3);
CommonFeasibleRate=100*mean(CommonFeasible,1).';
FeaRate=squeeze(100*mean(Feasible,1));
[Pc_dBm,Pc_lo,Pc_hi,Pc_n,Pc_mean_W]= ...
    masked_mean_ci_to_db(Pc,CommonFeasible,1e-3);
[Theta_mu,Theta_lo,Theta_hi,Theta_n]= ...
    masked_mean_ci_linear(Theta,CommonFeasible);

P1=Pc(:,:,1); P2=Pc(:,:,2); P3=Pc(:,:,3); P4=Pc(:,:,4);
[GainRISFixed,GainRISFixedLo,GainRISFixedHi,GainRISFixedN]= ...
    paired_db_statistics(P1,P2,Feasible(:,:,1)&Feasible(:,:,2));
[GainRISRET,GainRISRETLo,GainRISRETHi,GainRISRETN]= ...
    paired_db_statistics(P3,P4,Feasible(:,:,3)&Feasible(:,:,4));
[GainRET,GainRETLo,GainRETHi,GainRETN]= ...
    paired_db_statistics(P1,P3,Feasible(:,:,1)&Feasible(:,:,3));
[GainJoint,GainJointLo,GainJointHi,GainJointN]= ...
    paired_db_statistics(P1,P4,Feasible(:,:,1)&Feasible(:,:,4));
[Synergy,SynergyLo,SynergyHi,SynergyN]=paired_db_statistics( ...
    P2.*P3,P1.*P4,CommonFeasible);

if SAVE_MAIN_FIGURE
%% Figure 1: fair common-feasible power comparison
fig1=figure('Position',[60 60 980 650],'Color','w'); hold on;
for s=1:S
    errorbar(gamma_dB_range,Pc_dBm(:,s),Pc_dBm(:,s)-Pc_lo(:,s), ...
        Pc_hi(:,s)-Pc_dBm(:,s),markers{s},'LineWidth',1.5,'MarkerSize',8);
end
grid on; xlabel('SINR threshold \gamma (dB)');
ylabel('Minimum communication power P_c (dBm)');
title(sprintf('Communication Power on Common-Feasible Realizations (minimum n = %d)', ...
    min(Pc_n(:)))); legend(labels,'Location','northwest');
if SAVE_MAIN_FIGURE; save_exp1_figure(fig1,figure_dir,'01_Minimum_Communication_Power_Common_Feasible'); end

end

if SAVE_SUPP_FIGURES
%% Figure 2: RET angle with hardware reference
fig2=figure('Position',[80 80 1120 480],'Color','w');
idxs=[3 4]; names={'No RIS','Optimized RIS'};
for jj=1:2
    s=idxs(jj); subplot(1,2,jj); hold on;
    errorbar(gamma_dB_range,Theta_mu(:,s),Theta_mu(:,s)-Theta_lo(:,s), ...
        Theta_hi(:,s)-Theta_mu(:,s),'-s','LineWidth',1.5,'MarkerSize',7);
    yline(cfg.theta_ref,'--k',sprintf('Fixed baseline %.0f^\circ',cfg.theta_ref));
    grid on; xlabel('\gamma (dB)'); ylabel('Optimal RET tilt (deg)'); title(names{jj});
end
sgtitle(sprintf('Monte Carlo Mean of %.0f^\circ-Grid RET Decisions',cfg.ret_angle_resolution_deg));
if SAVE_SUPP_FIGURES; save_exp1_figure(fig2,figure_dir,'02_Optimal_RET_Tilt_Angles_With_CI'); end

%% Figure 3: paired marginal gains
fig3=figure('Position',[100 100 900 590],'Color','w'); hold on;
errorbar(gamma_dB_range,GainRISFixed,GainRISFixed-GainRISFixedLo, ...
    GainRISFixedHi-GainRISFixed,'-o','LineWidth',1.6,'MarkerSize',8);
errorbar(gamma_dB_range,GainRISRET,GainRISRET-GainRISRETLo, ...
    GainRISRETHi-GainRISRET,'-s','LineWidth',1.6,'MarkerSize',8);
errorbar(gamma_dB_range,GainRET,GainRET-GainRETLo, ...
    GainRETHi-GainRET,'-d','LineWidth',1.5,'MarkerSize',7);
grid on; xlabel('SINR threshold \gamma (dB)'); ylabel('Paired power saving (dB)');
legend({'RIS gain at fixed RET','RIS gain after RET optimization','RET-only gain'}, ...
    'Location','best'); title('Realization-Paired Marginal Communication-Power Gains');
if SAVE_SUPP_FIGURES; save_exp1_figure(fig3,figure_dir,'03_Paired_Marginal_Power_Saving'); end

%% Figure 4: additivity/synergy diagnostic
fig4=figure('Position',[120 120 900 560],'Color','w');
errorbar(gamma_dB_range,Synergy,Synergy-SynergyLo,SynergyHi-Synergy, ...
    '-p','LineWidth',1.7,'MarkerSize',9); hold on; yline(0,'--k'); grid on;
xlabel('SINR threshold \gamma (dB)');
ylabel('\Delta_{syn}=G_{joint}-G_{RIS}-G_{RET} (dB)');
title('RET-RIS Super-Additivity Diagnostic (Positive Means Super-Additive)');
if SAVE_SUPP_FIGURES; save_exp1_figure(fig4,figure_dir,'04_RET_RIS_Synergy_Index'); end

end

%% Tables and exact data
vnames=[{'SINR_dB'},arrayfun(@(x)sprintf('Mean_%d',x),1:S,'UniformOutput',false), ...
    arrayfun(@(x)sprintf('CI95Low_%d',x),1:S,'UniformOutput',false), ...
    arrayfun(@(x)sprintf('CI95High_%d',x),1:S,'UniformOutput',false), ...
    arrayfun(@(x)sprintf('Ncommon_%d',x),1:S,'UniformOutput',false), ...
    {'CommonFeasible_percent'}, ...
    arrayfun(@(x)sprintf('FeasiblePercent_%d',x),1:S,'UniformOutput',false)];
Tpower=array2table([gamma_dB_range(:),Pc_dBm,Pc_lo,Pc_hi,Pc_n, ...
    CommonFeasibleRate,FeaRate], 'VariableNames',vnames);
Tgain=table(gamma_dB_range(:),GainRISFixed,GainRISFixedLo,GainRISFixedHi, ...
    GainRISFixedN,GainRISRET,GainRISRETLo,GainRISRETHi,GainRISRETN, ...
    GainRET,GainRETLo,GainRETHi,GainRETN,GainJoint,GainJointLo,GainJointHi, ...
    GainJointN,Synergy,SynergyLo,SynergyHi,SynergyN, ...
    'VariableNames',{'SINR_dB','RIS_FixedRET_Mean_dB','RIS_FixedRET_Low_dB', ...
    'RIS_FixedRET_High_dB','RIS_FixedRET_N','RIS_OptimizedRET_Mean_dB', ...
    'RIS_OptimizedRET_Low_dB','RIS_OptimizedRET_High_dB','RIS_OptimizedRET_N', ...
    'RET_Only_Mean_dB','RET_Only_Low_dB','RET_Only_High_dB','RET_Only_N', ...
    'Joint_Mean_dB','Joint_Low_dB','Joint_High_dB','Joint_N', ...
    'Synergy_Mean_dB','Synergy_Low_dB','Synergy_High_dB','Synergy_N'});
disp(Tpower); disp(Tgain);
results_dir=fullfile(root,'results'); if ~exist(results_dir,'dir'); mkdir(results_dir); end
matfile=fullfile(results_dir,sprintf('exp1_revised_%s.mat',lower(RUN_MODE)));
save(matfile,'cfg','labels','gamma_dB_range','Pc','Feasible','Theta', ...
    'CommonFeasible','CommonFeasibleRate','FeaRate','Pc_dBm','Pc_lo','Pc_hi', ...
    'Pc_n','Pc_mean_W','Theta_mu','Theta_lo','Theta_hi','Theta_n', ...
    'GainRISFixed','GainRISFixedLo','GainRISFixedHi','GainRISRET', ...
    'GainRISRETLo','GainRISRETHi','GainRET','GainRETLo','GainRETHi', ...
    'GainJoint','GainJointLo','GainJointHi','Synergy','SynergyLo','SynergyHi', ...
    'Tpower','Tgain','runtime');
writetable(Tpower,fullfile(results_dir,sprintf('exp1_power_statistics_%s.csv',lower(RUN_MODE))));
writetable(Tgain,fullfile(results_dir,sprintf('exp1_paired_gains_%s.csv',lower(RUN_MODE))));
fprintf('\nExp1 completed in %.1f s. Common-feasible minimum count: %d/%d.\n', ...
    runtime,min(Pc_n(:)),ITER);
if SAVE_FIGURES; fprintf('Figures: %s\n',figure_dir); else; fprintf('Per-experiment figures: not exported (%s profile)\n',FIG_POLICY.profile); end
fprintf('MAT: %s\n',matfile);
fprintf('Total script runtime: %.1f s\n',toc(t_total));

function save_exp1_figure(fig,folder,name)
try
    exportgraphics(fig,fullfile(folder,[name '.png']),'Resolution',300);
catch
    print(fig,fullfile(folder,[name '.png']),'-dpng','-r300');
end
savefig(fig,fullfile(folder,[name '.fig']));
end

function exp1_progress(total)
persistent count t0
if isempty(count)||count>=total; count=0; t0=tic; end
count=count+1; elapsed=toc(t0); eta=elapsed/max(count,1)*(total-count);
fprintf('\rExp1 progress: %3.0f%% (%d/%d), elapsed %.1f min, ETA %.1f min', ...
    100*count/total,count,total,elapsed/60,eta/60);
if count>=total; fprintf('\n'); end
end
