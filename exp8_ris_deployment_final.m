% =========================================================================
% Experiment 8 (): RIS Deployment Position Sensitivity
% Compares no-RIS, fixed-RET/RIS radar-oriented, and joint RET-RIS designs
% on paired roadside realizations while moving the RIS along the road.
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
root=fileparts(mfilename('fullpath')); if isempty(root); root=pwd; end
addpath(root); addpath(fullfile(root,'function'));
if exist('cvx_begin','file')~=2; error('CVX is required. Run cvx_setup first.'); end

ITER=cfg.ITER_EXP8; xRIS=cfg.ris_x_range_exp8(:).'; P=numel(xRIS);
gamma=10^(cfg.gamma_exp8_dB/10); P_total=cfg.P_total_exp8;
labels={'Fixed RET, No RIS','Fixed RET, RIS-RadarSNR','Proposed RET-RIS-RadarSNR'};
S=numel(labels);
Obj=nan(ITER,P,S); Pcomm=nan(ITER,P,S); Feasible=false(ITER,P,S);
Theta=nan(ITER,P,S); Runtime=nan(ITER,P,S);

fprintf('\nExperiment 8  RIS deployment position (%s mode, MC=%d)\n',upper(RUN_MODE),ITER);
use_parallel=setup_parallel(cfg,ITER,root);
if use_parallel
    parfor mc=1:ITER
        [o,p,f,th,rt]=run_one(mc,cfg,xRIS,gamma,P_total);
        Obj(mc,:,:)=o; Pcomm(mc,:,:)=p; Feasible(mc,:,:)=f;
        Theta(mc,:,:)=th; Runtime(mc,:,:)=rt;
    end
else
    for mc=1:ITER
        [o,p,f,th,rt]=run_one(mc,cfg,xRIS,gamma,P_total);
        Obj(mc,:,:)=o; Pcomm(mc,:,:)=p; Feasible(mc,:,:)=f;
        Theta(mc,:,:)=th; Runtime(mc,:,:)=rt;
        if SHOW_PROGRESS; fprintf('Exp8 %d/%d\n',mc,ITER); end
    end
end

Common=all(Feasible,3);
[muObj,loObj,hiObj,nObj]=masked_mean_ci_to_db(Obj,Common,1);
[muPc,loPc,hiPc]=masked_mean_ci_to_db(Pcomm,Common,1e-3);
FeaRate=100*squeeze(mean(Feasible,1));
Gain=nan(ITER,P,2);
for pp=1:P
    ok=Common(:,pp)&Obj(:,pp,1)>0&Obj(:,pp,2)>0&Obj(:,pp,3)>0;
    Gain(ok,pp,1)=10*log10(Obj(ok,pp,3)./Obj(ok,pp,1));
    Gain(ok,pp,2)=10*log10(Obj(ok,pp,3)./Obj(ok,pp,2));
end
[muGain,loGain,hiGain]=masked_mean_ci_linear(Gain,isfinite(Gain));
[thetaMu,thetaLo,thetaHi]=masked_mean_ci_linear(Theta,Feasible);

if SAVE_MAIN_FIGURE
run_stamp=datestr(now,'yyyymmdd_HHMMSS');
figdir=fullfile(root,sprintf('Exp8_Figures_%s_%s',lower(RUN_MODE),run_stamp));
if SAVE_FIGURES && ~exist(figdir,'dir'); mkdir(figdir); end
fig=figure('Color','w','Position',[50 50 1450 900]); tl=tiledlayout(2,2,'TileSpacing','compact');
nexttile; hold on;
for s=1:S; errorbar(xRIS,muObj(:,s),muObj(:,s)-loObj(:,s),hiObj(:,s)-muObj(:,s),'-o','LineWidth',1.6); end
grid on; xlabel('RIS longitudinal position x_{RIS} (m)'); ylabel('Weighted radar SNR (dB)'); legend(labels,'Location','best'); title('(a) Deployment-dependent sensing performance');
nexttile; hold on;
for s=1:S; errorbar(xRIS,muPc(:,s),muPc(:,s)-loPc(:,s),hiPc(:,s)-muPc(:,s),'-s','LineWidth',1.5); end
grid on; xlabel('RIS longitudinal position x_{RIS} (m)'); ylabel('Communication power (dBm)'); title('(b) Communication cost on common samples');
nexttile; hold on;
errorbar(xRIS,muGain(:,1),muGain(:,1)-loGain(:,1),hiGain(:,1)-muGain(:,1),'-o','LineWidth',1.6);
errorbar(xRIS,muGain(:,2),muGain(:,2)-loGain(:,2),hiGain(:,2)-muGain(:,2),'-s','LineWidth',1.6);
yline(0,'--k'); grid on; xlabel('RIS longitudinal position x_{RIS} (m)'); ylabel('Proposed radar-SNR gain (dB)'); legend({'over no RIS','over fixed-RET RIS'},'Location','best'); title('(c) Paired deployment gain');
nexttile; hold on;
for s=2:S; errorbar(xRIS,thetaMu(:,s),thetaMu(:,s)-thetaLo(:,s),thetaHi(:,s)-thetaMu(:,s),'-d','LineWidth',1.5); end
grid on; xlabel('RIS longitudinal position x_{RIS} (m)'); ylabel('Selected RET tilt (deg)'); legend(labels(2:3),'Location','best'); title('(d) RET adaptation to deployment');
sgtitle(tl,sprintf('RIS Deployment Sensitivity, QoS=%.1f dB, P_T=%.1f dBm',cfg.gamma_exp8_dB,cfg.P_total_exp8_dBm));
if SAVE_MAIN_FIGURE; save_figure(fig,figdir,'01_RIS_Deployment_Sensitivity'); end

end

results_dir=fullfile(root,'results'); if ~exist(results_dir,'dir'); mkdir(results_dir); end
T=table(); rowCount=0;
for pp=1:P
    for s=1:S
        row=table(xRIS(pp),string(labels{s}),muObj(pp,s),loObj(pp,s),hiObj(pp,s), ...
            muPc(pp,s),FeaRate(pp,s),nObj(pp,s),'VariableNames', ...
            {'RIS_x_m','Scheme','RadarSNR_Mean_dB','CI95_Low_dB','CI95_High_dB', ...
             'Pcomm_Mean_dBm','Feasibility_percent','Common_samples'});
        rowCount=rowCount+1;
        if rowCount==1; T=row; else; T=[T;row]; end %#ok<AGROW>
    end
end
matfile=fullfile(results_dir,sprintf('exp8_ris_deployment_%s.mat',lower(RUN_MODE)));
save(matfile,'cfg','xRIS','labels','Obj','Pcomm','Feasible','Theta','Runtime','Common', ...
    'muObj','loObj','hiObj','muPc','loPc','hiPc','Gain','muGain','loGain','hiGain','T');
writetable(T,fullfile(results_dir,sprintf('exp8_ris_deployment_%s.csv',lower(RUN_MODE))));
fprintf('Experiment 8 complete in %.1f s. Results: %s\n',toc(t_total),matfile);

function [Obj,Pc,Fea,Theta,Runtime]=run_one(mc,cfg,xRIS,gamma,P_total)
P=numel(xRIS); S=3; Obj=nan(1,P,S); Pc=Obj; Fea=false(1,P,S); Theta=Obj; Runtime=Obj;
seed=cfg.seed_base+cfg.seed_stride*8+mc;
rng(seed+19,'twister');
xu=cfg.road_x_min+(cfg.road_x_max-cfg.road_x_min)*rand(1,cfg.K);
yu=-cfg.road_y_half+2*cfg.road_y_half*rand(1,cfg.K);
for pp=1:P
    cp=cfg; cp.ris_x=xRIS(pp); cp.ris_y=cfg.ris_y_exp8; cp=refresh_ret_ris_geometry(cp);
    scen=generate_isac_scenario(cp,cp.N,cp.K,seed,xu,yu,struct());
    rng(seed+10000,'twister'); out=evaluate_ret_ris_schemes(cp,scen,gamma,P_total);
    Obj(1,pp,1)=out.Obj(1); Pc(1,pp,1)=out.Pcomm(1); Fea(1,pp,1)=out.Feasible(1); Theta(1,pp,1)=out.Theta(1); Runtime(1,pp,1)=out.RuntimeCumulative(1);
    Obj(1,pp,3)=out.Obj(5); Pc(1,pp,3)=out.Pcomm(5); Fea(1,pp,3)=out.Feasible(5); Theta(1,pp,3)=out.Theta(5); Runtime(1,pp,3)=out.RuntimeCumulative(5);
    cfix=cp; cfix.exp6_fresh_candidate_mode='full'; cfix.exp6_fresh_allow_random=true;
    rng(seed+20000,'twister'); tt=tic;
    [~,rr]=optimize_fixed_ret_ris_radar(cfix,scen,cp.theta_ref,ones(cp.N,1),gamma,P_total);
    Runtime(1,pp,2)=toc(tt);
    if rr.feasible; Obj(1,pp,2)=rr.obj_snr; Pc(1,pp,2)=rr.p_comm; Fea(1,pp,2)=true; Theta(1,pp,2)=rr.theta_deg; end
end
end

function use_parallel=setup_parallel(cfg,ITER,root)
use_parallel=false;
try
    pool=gcp('nocreate'); if isempty(pool); parpool(min(cfg.num_workers,ITER)); end
    pctRunOnAll(['addpath(''' root ''');']); pctRunOnAll(['addpath(''' fullfile(root,'function') ''');']);
    use_parallel=true;
catch ME
    warning('Parallel execution unavailable: %s. Serial fallback is used.',ME.message);
end
end

function save_figure(fig,folder,name)
try; exportgraphics(fig,fullfile(folder,[name '.png']),'Resolution',300); catch; print(fig,fullfile(folder,[name '.png']),'-dpng','-r300'); end
savefig(fig,fullfile(folder,[name '.fig']));
end
