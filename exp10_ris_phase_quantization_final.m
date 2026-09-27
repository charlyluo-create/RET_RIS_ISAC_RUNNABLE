% =========================================================================
% Experiment 10 (): Practical RIS Phase Quantization
% The continuous joint RET-RIS solution is quantized to 1/2/3-bit phases.
% Active beams are then re-optimized under the quantized hardware setting.
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

ITER=cfg.ITER_EXP10; bits=cfg.ris_phase_bits_exp10(:).'; B=numel(bits);
gamma=10^(cfg.gamma_exp10_dB/10); P_total=cfg.P_total_exp10;
labels=cell(1,B); for b=1:B; if isinf(bits(b)); labels{b}='Continuous'; else; labels{b}=sprintf('%d-bit',bits(b)); end; end
Obj=nan(ITER,B); Pcomm=nan(ITER,B); Feasible=false(ITER,B); PhaseRMSE_deg=nan(ITER,B);
NoRISObj=nan(ITER,1); NoRISFeasible=false(ITER,1); Theta=nan(ITER,1);

fprintf('\nExperiment 10  RIS phase quantization (%s mode, MC=%d)\n',upper(RUN_MODE),ITER);
use_parallel=setup_parallel(cfg,ITER,root);
if use_parallel
    parfor mc=1:ITER
        [Obj(mc,:),Pcomm(mc,:),Feasible(mc,:),PhaseRMSE_deg(mc,:), ...
            NoRISObj(mc),NoRISFeasible(mc),Theta(mc)]=run_one(mc,cfg,bits,gamma,P_total);
    end
else
    for mc=1:ITER
        [Obj(mc,:),Pcomm(mc,:),Feasible(mc,:),PhaseRMSE_deg(mc,:), ...
            NoRISObj(mc),NoRISFeasible(mc),Theta(mc)]=run_one(mc,cfg,bits,gamma,P_total);
        if SHOW_PROGRESS; fprintf('Exp10 %d/%d\n',mc,ITER); end
    end
end

Common=all(Feasible,2)&NoRISFeasible;
mask=repmat(Common,1,B);
[muObj,loObj,hiObj,nObj]=masked_mean_ci_to_db(Obj,mask,1);
[muPc,loPc,hiPc]=masked_mean_ci_to_db(Pcomm,mask,1e-3);
[noMu,noLo,noHi]=masked_mean_ci_to_db(NoRISObj,Common,1);
Loss=nan(ITER,B);
for b=1:B
    ok=Common&Obj(:,b)>0&Obj(:,end)>0;
    Loss(ok,b)=10*log10(Obj(ok,end)./Obj(ok,b));
end
[lossMu,lossLo,lossHi]=masked_mean_ci_linear(Loss,isfinite(Loss));
[phaseMu,phaseLo,phaseHi]=masked_mean_ci_linear(PhaseRMSE_deg,Feasible);
FeaRate=100*mean(Feasible,1);

x=1:B;
if SAVE_MAIN_FIGURE
run_stamp=datestr(now,'yyyymmdd_HHMMSS');
figdir=fullfile(root,sprintf('Exp10_Figures_%s_%s',lower(RUN_MODE),run_stamp));
if SAVE_FIGURES && ~exist(figdir,'dir'); mkdir(figdir); end
fig=figure('Color','w','Position',[70 70 1450 820]); tl=tiledlayout(2,2,'TileSpacing','compact');
nexttile;
errorbar(x,muObj,muObj-loObj,hiObj-muObj,'-o','LineWidth',1.7,'MarkerSize',8); hold on;
yline(noMu,'--','No-RIS reference'); grid on; set(gca,'XTick',x,'XTickLabel',labels); ylabel('Weighted radar SNR (dB)'); title('(a) Quantized-phase sensing performance');
nexttile;
errorbar(x,lossMu,lossMu-lossLo,lossHi-lossMu,'-s','LineWidth',1.7,'MarkerSize',8); yline(0,'--k'); grid on; set(gca,'XTick',x,'XTickLabel',labels); ylabel('Loss relative to continuous phase (dB)'); title('(b) Paired quantization loss');
nexttile;
errorbar(x,phaseMu,phaseMu-phaseLo,phaseHi-phaseMu,'-d','LineWidth',1.6,'MarkerSize',8); grid on; set(gca,'XTick',x,'XTickLabel',labels); ylabel('Phase RMSE (deg)'); title('(c) Hardware phase distortion');
nexttile; yyaxis left;
errorbar(x,muPc,muPc-loPc,hiPc-muPc,'-^','LineWidth',1.5); ylabel('Communication power (dBm)'); yyaxis right;
plot(x,FeaRate,'-o','LineWidth',1.5); ylabel('Feasibility (%)'); ylim([0 105]); grid on; set(gca,'XTick',x,'XTickLabel',labels); title('(d) Communication cost and feasibility');
sgtitle(tl,sprintf('RIS Phase Quantization, QoS=%.1f dB, P_T=%.1f dBm',cfg.gamma_exp10_dB,cfg.P_total_exp10_dBm));
if SAVE_MAIN_FIGURE; save_figure(fig,figdir,'01_RIS_Phase_Quantization'); end

end

results_dir=fullfile(root,'results'); if ~exist(results_dir,'dir'); mkdir(results_dir); end
T=table(bits(:),string(labels(:)),muObj(:),loObj(:),hiObj(:),lossMu(:), ...
    phaseMu(:),muPc(:),FeaRate(:),nObj(:),'VariableNames', ...
    {'Phase_bits','Label','RadarSNR_Mean_dB','CI95_Low_dB','CI95_High_dB', ...
     'Loss_vs_continuous_dB','Phase_RMSE_deg','Pcomm_Mean_dBm','Feasibility_percent','Common_samples'});
matfile=fullfile(results_dir,sprintf('exp10_ris_quantization_%s.mat',lower(RUN_MODE)));
save(matfile,'cfg','bits','labels','Obj','Pcomm','Feasible','PhaseRMSE_deg','NoRISObj', ...
    'NoRISFeasible','Theta','Common','muObj','loObj','hiObj','Loss','lossMu','lossLo','lossHi','T');
writetable(T,fullfile(results_dir,sprintf('exp10_ris_quantization_%s.csv',lower(RUN_MODE))));
fprintf('Experiment 10 complete in %.1f s. Results: %s\n',toc(t_total),matfile);

function [Obj,Pc,Fea,RMSE,noObj,noFea,theta]=run_one(mc,cfg,bits,gamma,P_total)
B=numel(bits); Obj=nan(1,B); Pc=Obj; Fea=false(1,B); RMSE=Obj;
noObj=NaN; noFea=false; theta=NaN;
seed=cfg.seed_base+cfg.seed_stride*10+mc;
scen=generate_isac_scenario(cfg,cfg.N,cfg.K,seed,[],[],struct());
rng(seed+100000,'twister'); out=evaluate_ret_ris_schemes(cfg,scen,gamma,P_total);
if out.Feasible(1); noObj=out.Obj(1); noFea=true; end
if ~out.Feasible(5); return; end
theta=out.Theta(5); phi=out.Phi{5};
for b=1:B
    phiq=quantize_ris_phase(phi,bits(b));
    rr=evaluate_fixed_ret_ris_config(cfg,scen,theta,phiq,gamma,P_total,true);
    if rr.feasible
        Obj(b)=rr.obj_snr; Pc(b)=rr.p_comm; Fea(b)=true;
        d=angle(phiq.*conj(phi)); RMSE(b)=sqrt(mean(d.^2))*180/pi;
    end
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
