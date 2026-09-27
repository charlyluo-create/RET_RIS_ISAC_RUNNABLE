% =========================================================================
% Experiment 11 (): RET Angle Resolution versus Performance and Runtime
% A fixed physical local-search span is maintained across hardware grids.
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

ITER=cfg.ITER_EXP11; resRange=cfg.ret_resolution_deg_exp11(:).'; R=numel(resRange);
gamma=10^(cfg.gamma_exp11_dB/10); P_total=cfg.P_total_exp11;
Obj=nan(ITER,R); Pcomm=nan(ITER,R); Feasible=false(ITER,R); Runtime=nan(ITER,R); Theta=nan(ITER,R);
GridSize=nan(1,R); LocalCandidates=nan(1,R);
for r=1:R
    GridSize(r)=numel(0:resRange(r):35);
    LocalCandidates(r)=2*max(1,round(cfg.ret_local_search_span_deg_exp11/resRange(r)))+1;
end

fprintf('\nExperiment 11  RET angle resolution (%s mode, MC=%d)\n',upper(RUN_MODE),ITER);
use_parallel=setup_parallel(cfg,ITER,root);
if use_parallel
    parfor mc=1:ITER
        [Obj(mc,:),Pcomm(mc,:),Feasible(mc,:),Runtime(mc,:),Theta(mc,:)]=run_one(mc,cfg,resRange,gamma,P_total);
    end
else
    for mc=1:ITER
        [Obj(mc,:),Pcomm(mc,:),Feasible(mc,:),Runtime(mc,:),Theta(mc,:)]=run_one(mc,cfg,resRange,gamma,P_total);
        if SHOW_PROGRESS; fprintf('Exp11 %d/%d\n',mc,ITER); end
    end
end

Common=all(Feasible,2); mask=repmat(Common,1,R);
[muObj,loObj,hiObj,nObj]=masked_mean_ci_to_db(Obj,mask,1);
[muPc,loPc,hiPc]=masked_mean_ci_to_db(Pcomm,mask,1e-3);
[thetaMu,thetaLo,thetaHi]=masked_mean_ci_linear(Theta,mask);
[runMed,runQ1,runQ3]=median_iqr(Runtime);
Loss=nan(ITER,R);
for r=1:R
    ok=Common&Obj(:,r)>0&Obj(:,end)>0;
    Loss(ok,r)=10*log10(Obj(ok,end)./Obj(ok,r));
end
[lossMu,lossLo,lossHi]=masked_mean_ci_linear(Loss,isfinite(Loss));

x=1:R; labels=arrayfun(@(v)sprintf('%.1f°',v),resRange,'UniformOutput',false);
if SAVE_MAIN_FIGURE
run_stamp=datestr(now,'yyyymmdd_HHMMSS');
figdir=fullfile(root,sprintf('Exp11_Figures_%s_%s',lower(RUN_MODE),run_stamp));
if SAVE_FIGURES && ~exist(figdir,'dir'); mkdir(figdir); end
fig=figure('Color','w','Position',[70 70 1450 820]); tl=tiledlayout(2,2,'TileSpacing','compact');
nexttile; errorbar(x,muObj,muObj-loObj,hiObj-muObj,'-o','LineWidth',1.7,'MarkerSize',8); grid on; set(gca,'XTick',x,'XTickLabel',labels); ylabel('Weighted radar SNR (dB)'); title('(a) RET resolution performance');
nexttile; errorbar(x,lossMu,lossMu-lossLo,lossHi-lossMu,'-s','LineWidth',1.7,'MarkerSize',8); yline(0,'--k'); grid on; set(gca,'XTick',x,'XTickLabel',labels); ylabel('Loss relative to 0.5° grid (dB)'); title('(b) Paired resolution loss');
nexttile; errorbar(x,runMed,runMed-runQ1,runQ3-runMed,'-d','LineWidth',1.7,'MarkerSize',8); grid on; set(gca,'XTick',x,'XTickLabel',labels); ylabel('Complete runtime (s)'); title('(c) Performance-complexity cost');
nexttile; errorbar(x,thetaMu,thetaMu-thetaLo,thetaHi-thetaMu,'-^','LineWidth',1.6,'MarkerSize',8); grid on; set(gca,'XTick',x,'XTickLabel',labels); ylabel('Selected RET tilt (deg)'); title('(d) Selected hardware setting');
sgtitle(tl,sprintf('RET Angle Resolution, QoS=%.1f dB, P_T=%.1f dBm',cfg.gamma_exp11_dB,cfg.P_total_exp11_dBm));
if SAVE_MAIN_FIGURE; save_figure(fig,figdir,'01_RET_Resolution_Tradeoff'); end

end

results_dir=fullfile(root,'results'); if ~exist(results_dir,'dir'); mkdir(results_dir); end
T=table(resRange(:),GridSize(:),LocalCandidates(:),muObj(:),loObj(:),hiObj(:), ...
    lossMu(:),runMed(:),runQ1(:),runQ3(:),thetaMu(:),muPc(:),nObj(:), ...
    'VariableNames',{'RET_resolution_deg','Full_grid_points','Local_candidates', ...
    'RadarSNR_Mean_dB','CI95_Low_dB','CI95_High_dB','Loss_vs_0p5deg_dB', ...
    'Runtime_median_s','Runtime_Q1_s','Runtime_Q3_s','Selected_tilt_mean_deg', ...
    'Pcomm_Mean_dBm','Common_samples'});
matfile=fullfile(results_dir,sprintf('exp11_ret_resolution_%s.mat',lower(RUN_MODE)));
save(matfile,'cfg','resRange','Obj','Pcomm','Feasible','Runtime','Theta','GridSize', ...
    'LocalCandidates','Common','muObj','loObj','hiObj','Loss','lossMu','lossLo','lossHi','T');
writetable(T,fullfile(results_dir,sprintf('exp11_ret_resolution_%s.csv',lower(RUN_MODE))));
fprintf('Experiment 11 complete in %.1f s. Results: %s\n',toc(t_total),matfile);

function [Obj,Pc,Fea,Run,Theta]=run_one(mc,cfg,resRange,gamma,P_total)
R=numel(resRange); Obj=nan(1,R); Pc=Obj; Fea=false(1,R); Run=Obj; Theta=Obj;
seed=cfg.seed_base+cfg.seed_stride*11+mc;
rng(seed+31,'twister'); xu=cfg.road_x_min+(cfg.road_x_max-cfg.road_x_min)*rand(1,cfg.K); yu=-cfg.road_y_half+2*cfg.road_y_half*rand(1,cfg.K);
for r=1:R
    cr=cfg; cr.ret_angle_resolution_deg=resRange(r); cr.theta_tilt_range=0:resRange(r):35;
    cr.tilt_local_half_window=max(1,round(cr.ret_local_search_span_deg_exp11/resRange(r)));
    cr=refresh_ret_ris_geometry(cr);
    scen=generate_isac_scenario(cr,cr.N,cr.K,seed,xu,yu,struct());
    rng(seed+110000,'twister'); tt=tic; out=evaluate_ret_ris_schemes(cr,scen,gamma,P_total); Run(r)=toc(tt);
    if out.Feasible(5); Obj(r)=out.Obj(5); Pc(r)=out.Pcomm(5); Fea(r)=true; Theta(r)=out.Theta(5); end
end
end

function [med,q1,q3]=median_iqr(X)
med=nan(1,size(X,2)); q1=med; q3=med;
for j=1:size(X,2); v=X(:,j); v=v(isfinite(v)); if ~isempty(v); med(j)=median(v); q1(j)=empirical_quantile(v,0.25); q3(j)=empirical_quantile(v,0.75); end; end
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
