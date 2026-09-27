function output_dir = generate_paper_figures(run_mode)
%GENERATE_PAPER_FIGURES Build six composite manuscript figures with consistent formatting.
% Updated to use slanted categorical x-axis labels for readability.
%
% The 12 experiment scripts remain available as data/diagnostic generators,
% but the manuscript uses only six experimental composite figures:
%   F1  Communication-sensing coupling       (Exp. 1--3)
%   F2  Mechanism and same-objective ablation (Exp. 4)
%   F3  Scalability and RIS deployment        (Exp. 5 and 8)
%   F4  CSI sensitivity and hardware limits    (Exp. 9--11)
%   F5  Mobility and RIS staleness             (Exp. 6)
%   F6  Algorithm validation                   (Exp. 7; Exp. 12 is tables)
%
% Output formats: vector PDF, 600-dpi PNG and editable FIG.
%
% Manuscript legend/style policy:
%   Fixed RET/no RIS     : blue, solid, circle
%   Joint/Fixed RET-RIS  : orange, dashed, square
%   Proposed MA-RCG-MF   : yellow, dash-dot, diamond
% The same scheme keeps the same visual identity across manuscript figures,
% and every multi-scheme panel includes an explicit legend so the figure is
% self-contained even in grayscale or when panels are viewed separately.

if nargin<1 || isempty(run_mode)
    env_mode=getenv('RET_RIS_RUN_MODE');
    if isempty(env_mode); run_mode='paper'; else; run_mode=lower(env_mode); end
end
run_mode=lower(char(run_mode));
root=fileparts(mfilename('fullpath')); if isempty(root); root=pwd; end
addpath(root); addpath(fullfile(root,'function'));
results_dir=fullfile(root,'results');
if ~exist(results_dir,'dir'); error('Results folder does not exist: %s',results_dir); end
stamp=datestr(now,'yyyymmdd_HHMMSS');
output_dir=fullfile(root,sprintf('Main_Figures_%s_%s',run_mode,stamp));
mkdir(output_dir);

D1=load_required(results_dir,sprintf('exp1_revised_%s.mat',run_mode));
D2=load_required(results_dir,sprintf('exp2_revised_%s.mat',run_mode));
D3=load_required(results_dir,sprintf('exp3_qos_revised_%s.mat',run_mode));
D4=load_required(results_dir,sprintf('exp4_ablation_revised_%s.mat',run_mode));
D5=load_required(results_dir,sprintf('exp5_scalability_revised_%s.mat',run_mode));
D6=load_required(results_dir,sprintf('exp6_multitimescale_mobility_%s.mat',run_mode));
cfgExpected=unified_ret_ris_isac_config(run_mode);
if ~isfield(D6,'cfg') || ~isfield(D6.cfg,'exp6_pairing_revision') || ...
        ~strcmp(D6.cfg.exp6_pairing_revision,cfgExpected.exp6_pairing_revision)
    error('F5 requires the strict-paired Exp.6 revision. Re-run the paired Exp.6 before building figures.');
end
D7=load_required(results_dir,sprintf('exp7_algorithm_convergence_%s.mat',run_mode));
D8=load_required(results_dir,sprintf('exp8_ris_deployment_%s.mat',run_mode));
D9=load_required(results_dir,sprintf('exp9_csi_uncertainty_%s.mat',run_mode));
D10=load_required(results_dir,sprintf('exp10_ris_quantization_%s.mat',run_mode));
D11=load_required(results_dir,sprintf('exp11_ret_resolution_%s.mat',run_mode));

%% Figure 1: communication-sensing coupling (Exp. 1--3)
fig=figure('Color','w','Position',[80 60 1400 1180]);
tl=tiledlayout(2,2,'TileSpacing','loose','Padding','loose');

nexttile; hold on;
idx1=1:min(4,size(D1.Pc_dBm,2));
for s=idx1
    h=errorbar(D1.gamma_dB_range,D1.Pc_dBm(:,s), ...
        D1.Pc_dBm(:,s)-D1.Pc_lo(:,s),D1.Pc_hi(:,s)-D1.Pc_dBm(:,s), ...
        'LineWidth',1.25,'MarkerSize',5);
    lab=compact_labels(D1.labels(s)); apply_label_style(h,lab{1});
end
grid on; xlabel('Communication SINR threshold (dB)');
ylabel('Minimum communication power (dBm)');
title('(a) Communication power');
legend(compact_labels(D1.labels(idx1)),'Location','northwest','FontSize',12.5,'Box','off');

nexttile; hold on;
idx2=[1 4 5]; idx2=idx2(idx2<=size(D2.Obj_dB,2));
for jj=1:numel(idx2)
    s=idx2(jj);
    h=errorbar(D2.range_P_dBm,D2.Obj_dB(:,s), ...
        D2.Obj_dB(:,s)-D2.Obj_lo_dB(:,s),D2.Obj_hi_dB(:,s)-D2.Obj_dB(:,s), ...
        'LineWidth',1.25,'MarkerSize',5);
    apply_three_scheme_style(h,jj);
end
grid on; xlabel('Total transmit power (dBm)'); ylabel('Weighted radar SNR (dB)');
title('(b) Sensing versus total power');
legend({'Fixed RET/no RIS','Joint RET-RIS MinPc','MA-RCG-MF'}, ...
    'Location','northwest','FontSize',12.5,'Box','off');

nexttile; hold on;
idx3=[1 4 5]; idx3=idx3(idx3<=size(D3.Obj_dB,2));
for jj=1:numel(idx3)
    s=idx3(jj);
    h=errorbar(D3.gamma_dB_range,D3.Obj_dB(:,s), ...
        D3.Obj_dB(:,s)-D3.Obj_lo_dB(:,s),D3.Obj_hi_dB(:,s)-D3.Obj_dB(:,s), ...
        'LineWidth',1.25,'MarkerSize',5);
    apply_three_scheme_style(h,jj);
end
grid on; xlabel('Communication SINR threshold (dB)'); ylabel('Weighted radar SNR (dB)');
title('(c) Sensing versus communication QoS');
legend({'Fixed RET/no RIS','Joint RET-RIS MinPc','MA-RCG-MF'}, ...
    'Location','southwest','FontSize',12.5,'Box','off');

nexttile; hold on;
rate3=zeros(numel(D3.gamma_dB_range),numel(idx3));
for jj=1:numel(idx3)
    s=idx3(jj); rate3(:,jj)=100*D3.Fea_rate(:,s);
    h=plot(D3.gamma_dB_range,rate3(:,jj),'LineWidth',1.25,'MarkerSize',5);
    apply_three_scheme_style(h,jj);
end
grid on; set_feasibility_ylim(gca,rate3); xlabel('Communication SINR threshold (dB)');
ylabel('Feasibility rate (%)'); title('(d) QoS feasibility');
legend({'Fixed RET/no RIS','Joint RET-RIS MinPc','MA-RCG-MF'}, ...
    'Location','southwest','FontSize',12.5,'Box','off');
% Joint RET-RIS MinPc and MA-RCG-MF can overlap exactly in this panel.
% Keep the axes uncluttered; explain the overlap in the manuscript caption.
sgtitle(tl,'Communication-Sensing Coupling','FontSize',18,'FontWeight','bold');
save_master(fig,output_dir,'F1_Communication_Sensing_Coupling'); close(fig);

%% Figure 2: dual mechanism and same-objective hardware ablation (Exp. 4)
S4=size(D4.Obj,2); mask4=repmat(D4.MechanismCommon,1,S4);
[pc4,pcl4,pch4]=masked_mean_ci_to_db(D4.Pcomm,mask4,1e-3);
[ge4,gel4,geh4]=masked_mean_ci_to_db(D4.G_eff,mask4,1);
[ob4,obl4,obh4]=masked_mean_ci_to_db(D4.Obj,mask4,1);
[pr4,prl4,prh4]=masked_mean_ci_linear(100*D4.Pradar/D4.P_total,mask4);
short4={'No RIS','Fixed/RIS','RET only','Joint','MA-RCG-MF'};
Ssame=size(D4.ObjSame,2); masksame=repmat(D4.SameCommon,1,Ssame);
[obs,obsl,obsh]=masked_mean_ci_to_db(D4.ObjSame,masksame,1);
same_labels={'No RIS','RET only','RIS only','Joint'};
fig=figure('Color','w','Position',[80 60 1400 1180]);
tl=tiledlayout(2,2,'TileSpacing','loose','Padding','loose');
nexttile; bar_ci(pr4,prl4,prh4,short4); ylabel('Residual sensing-power share (%)');
title('(a) Available sensing-power share'); grid on;
nexttile; bar_ci(ge4,gel4,geh4,short4); ylabel('Unit-power sensing gain (dB/W)');
title('(b) Per-watt sensing-channel gain'); grid on;
nexttile; bar_ci(ob4,obl4,obh4,short4); ylabel('Weighted radar SNR (dB)');
title('(c) Final sensing outcome'); grid on;
nexttile; bar_ci(obs,obsl,obsh,same_labels); ylabel('Weighted radar SNR (dB)');
title('(d) Same-objective hardware ablation'); grid on;
sgtitle(tl,'Mechanism and Hardware Ablation','FontSize',18,'FontWeight','bold');
save_master(fig,output_dir,'F2_Mechanism_and_Hardware_Ablation'); close(fig);

%% Figure 3: scalability and deployment (Exp. 5 and 8)
[objN,objNl,objNh]=masked_mean_ci_to_db(D5.ObjN,D5.CommonN,1);
[objK,objKl,objKh]=masked_mean_ci_to_db(D5.ObjK,D5.CommonK,1);
labels5={'Fixed RET/no RIS','Joint RET-RIS MinPc','MA-RCG-MF'};
labels8={'Fixed RET/no RIS','Fixed RET/RIS','MA-RCG-MF'};
fig=figure('Color','w','Position',[80 60 1400 1180]);
tl=tiledlayout(2,2,'TileSpacing','loose','Padding','loose');

% (a) Exp.5: all three curves use an explicit legend.
nexttile; hold on;
% No-RIS does not depend on N: draw it as an explicit reference line.
h=errorbar(D5.N_range,objN(:,1),objN(:,1)-objNl(:,1),objNh(:,1)-objN(:,1), ...
    'LineWidth',1.15,'MarkerSize',4); apply_reference_style(h);
for ss=2:size(objN,2)
    h=errorbar(D5.N_range,objN(:,ss),objN(:,ss)-objNl(:,ss),objNh(:,ss)-objN(:,ss), ...
        'LineWidth',1.25,'MarkerSize',5);
    apply_three_scheme_style(h,ss);
end
grid on; xlabel('RIS elements N'); ylabel('Radar SNR (dB)');
title('(a) RIS-size scalability');
legend({'No-RIS ref.','Joint MinPc','MA-RCG-MF'}, ...
    'Location','best','FontSize',14,'Box','off');

% (b) Exp.5: repeat the legend so this panel remains understandable when
% cropped or viewed independently.
nexttile; hold on;
for ss=1:size(objK,2)
    h=errorbar(D5.K_range,objK(:,ss),objK(:,ss)-objKl(:,ss),objKh(:,ss)-objK(:,ss), ...
        'LineWidth',1.25,'MarkerSize',5);
    apply_three_scheme_style(h,ss);
end
grid on; xlabel('Number of communication vehicles K'); ylabel('Radar SNR (dB)');
title('(b) Vehicular-load scalability');
legend(labels5,'Location','best','FontSize',14,'Box','off');

% (c) Exp.5: three feasibility curves.  Some curves can overlap exactly at
% 100%, therefore the legend is essential even when only two colors appear.
nexttile; hold on;
for ss=1:size(D5.FeaK_rate,2)
    h=plot(D5.K_range,D5.FeaK_rate(:,ss),'LineWidth',1.25,'MarkerSize',5);
    apply_three_scheme_style(h,ss);
end
grid on; set_feasibility_ylim(gca,D5.FeaK_rate); xlabel('Communication vehicles K'); ylabel('Feasibility rate (%)');
title('(c) High-load feasibility');
legend(labels5,'Location','southwest','FontSize',14,'Box','off');

% (d) Exp.8 uses a different second baseline from Exp.5.  Do NOT reuse the
% Exp.5 legend: the orange curve here is Fixed RET/RIS (radar-oriented), not
% Joint RET-RIS MinPc.
nexttile; hold on;
h=errorbar(D8.xRIS,D8.muObj(:,1),D8.muObj(:,1)-D8.loObj(:,1), ...
    D8.hiObj(:,1)-D8.muObj(:,1),'LineWidth',1.15,'MarkerSize',4); apply_reference_style(h);
for ss=2:size(D8.muObj,2)
    h=errorbar(D8.xRIS,D8.muObj(:,ss),D8.muObj(:,ss)-D8.loObj(:,ss), ...
        D8.hiObj(:,ss)-D8.muObj(:,ss),'LineWidth',1.25,'MarkerSize',5);
    apply_three_scheme_style(h,ss);
end
grid on; xlabel('RIS longitudinal position (m)'); ylabel('Radar SNR (dB)');
title('(d) Roadside RIS deployment');
legend({'No-RIS ref.','Fixed RET/RIS','MA-RCG-MF'}, ...
    'Location','best','FontSize',14,'Box','off');
sgtitle(tl,'Scalability and Roadside Deployment','FontSize',18,'FontWeight','bold');
save_master(fig,output_dir,'F3_Scalability_and_RIS_Deployment'); close(fig);

%% Figure 4: imperfect CSI and hardware resolution (Exp. 9--11)
fig=figure('Color','w','Position',[80 60 1400 1180]);
tl=tiledlayout(2,2,'TileSpacing','loose','Padding','loose');
assert(isfield(D9,'nmseTarget_dB'), ...
    'Exp.9 PAPER result is incompatible with the final clean release. Rerun Exp.9.');
x9=1:numel(D9.nmseTarget_dB);
x9labels=nmse_tick_labels(D9.nmseTarget_dB);
x9xlabel='Target CSI NMSE (dB)';
nexttile; hold on;
for ss=1:size(D9.muObj,2)
    h=errorbar(x9,D9.muObj(:,ss),D9.muObj(:,ss)-D9.loObj(:,ss), ...
        D9.hiObj(:,ss)-D9.muObj(:,ss),'LineWidth',1.25,'MarkerSize',5);
    apply_three_scheme_style(h,ss);
end
if ~isempty(x9labels); set(gca,'XTick',x9,'XTickLabel',x9labels,'XTickLabelRotation',25); xlim([0.7 numel(x9)+0.3]); end
grid on; xlabel(x9xlabel); ylabel('Outage-penalized radar SNR (dB)');
title('(a) Sensing under CSI mismatch');
legend({'Fixed RET/no RIS','Fixed RET/RIS','MA-RCG-MF'}, ...
    'Location','best','FontSize',14,'Box','off');
nexttile; hold on;
for ss=1:size(D9.QoSRate,2)
    h=plot(x9,D9.QoSRate(:,ss),'LineWidth',1.25,'MarkerSize',5);
    apply_three_scheme_style(h,ss);
end
if ~isempty(x9labels); set(gca,'XTick',x9,'XTickLabel',x9labels,'XTickLabelRotation',25); xlim([0.7 numel(x9)+0.3]); end
grid on; set_feasibility_ylim(gca,D9.QoSRate); xlabel(x9xlabel); ylabel('Actual QoS satisfaction rate (%)');
title('(b) Communication reliability under CSI mismatch');
% Perfect point is fixed by the cohort definition; explain this in the caption.
legend({'Fixed RET/no RIS','Fixed RET/RIS','MA-RCG-MF'}, ...
    'Location','southwest','FontSize',14,'Box','off');
nexttile;
x10=1:numel(D10.labels);
errorbar(x10,D10.lossMu(:),D10.lossMu(:)-D10.lossLo(:), ...
    D10.lossHi(:)-D10.lossMu(:),'-d','LineWidth',1.35,'MarkerSize',6);
grid on; yline(0,'--k','No loss','HandleVisibility','off');
set(gca,'XTick',x10,'XTickLabel',D10.labels); xlabel('RIS phase resolution');
ylabel('RIS phase loss (dB)'); title('(c) RIS phase quantization');
nexttile;
x11=1:numel(D11.resRange);
errorbar(x11,D11.lossMu(:),D11.lossMu(:)-D11.lossLo(:), ...
    D11.lossHi(:)-D11.lossMu(:),'-p','LineWidth',1.35,'MarkerSize',6);
ret_tick_labels=arrayfun(@(v) sprintf('%g^\\circ',v),D11.resRange,'UniformOutput',false);
grid on; yline(0,'--k','No loss','HandleVisibility','off');
set(gca,'XTick',x11,'XTickLabel',ret_tick_labels,'TickLabelInterpreter','tex');
xlabel('RET angle-grid resolution'); ylabel('RET-grid loss (dB)');
title('(d) RET resolution tradeoff');
sgtitle(tl,'CSI Sensitivity and Hardware Resolution','FontSize',18,'FontWeight','bold');
save_master(fig,output_dir,'F4_CSI_Sensitivity_and_Hardware'); close(fig);

%% Figure 5: multi-timescale mobility (Exp. 6)
% D6.update_ms is deliberately nonuniform.  Use equally spaced matrix columns
% and print the actual update intervals as tick labels; using imagesc(update_ms,
% ...) would visually interpolate nonuniform intervals as if they were uniform.
fig=figure('Color','w','Position',[80 100 1400 700]);
tl=tiledlayout(1,2,'TileSpacing','loose','Padding','loose');
climmax=max([D6.LossMean(:);D6.LossP90(:)],[],'omitnan');
if isempty(climmax)||~isfinite(climmax)||climmax<=0; climmax=1; end
x6=1:numel(D6.update_ms); x6labels=numeric_tick_labels(D6.update_ms);
y6=1:numel(D6.speed_kmh); y6labels=numeric_tick_labels(D6.speed_kmh);
nexttile; imagesc(x6,y6,D6.LossMean); axis xy; colorbar; caxis([0 climmax]);
set(gca,'XTick',x6,'XTickLabel',x6labels, ...
    'YTick',y6,'YTickLabel',y6labels);
xlim([0.5 numel(x6)+0.5]); ylim([0.5 numel(y6)+0.5]);
xlabel('RIS update period (ms)'); ylabel('Vehicle speed (km/h)');
title('(a) Mean paired staleness loss (dB)');
annotate_heatmap_values(gca,D6.LossMean,climmax);
nexttile; imagesc(x6,y6,D6.LossP90); axis xy; colorbar; caxis([0 climmax]);
set(gca,'XTick',x6,'XTickLabel',x6labels, ...
    'YTick',y6,'YTickLabel',y6labels);
xlim([0.5 numel(x6)+0.5]); ylim([0.5 numel(y6)+0.5]);
xlabel('RIS update period (ms)'); ylabel('Vehicle speed (km/h)');
title('(b) 90th-percentile paired loss (dB)');
annotate_heatmap_values(gca,D6.LossP90,climmax);
sgtitle(tl,'Strict-Paired Mobility and RIS Aging','FontSize',18,'FontWeight','bold');
save_master(fig,output_dir,'F5_Mobility_and_RIS_Staleness'); close(fig);

%% Figure 6: algorithm validation (Exp. 7; complexity/significance in tables)
idx=D7.main_idx(:).'; labels7=compact_labels(D7.labels(idx));
labels7Tick=algorithm_tick_labels(labels7);
[rtmed,rtq1,rtq3]=runtime_iqr(D7.Runtime);
fig=figure('Color','w','Position',[80 60 1400 1180]);
tl=tiledlayout(2,2,'TileSpacing','loose','Padding','loose');
nexttile; hold on;
for jj=1:numel(idx)
    s=idx(jj); yy=D7.ObjCommonMean_dB(s);
    h=errorbar(jj,yy,yy-D7.ObjCommonLo_dB(s),D7.ObjCommonHi_dB(s)-yy, ...
        'o','LineWidth',1.35,'MarkerSize',6);
    apply_label_style(h,labels7{jj});
end
grid on; xlim([0.35 numel(idx)+0.65]); set(gca,'XTick',1:numel(idx),'XTickLabel',labels7Tick,'XTickLabelRotation',28);
ylabel('Radar SNR (dB)'); title('(a) Same-model objective comparison');

nexttile; hold on;
runtimeVals=rtmed(idx); maLocal=find(idx==D7.proposed_idx,1);
for jj=1:numel(idx)
    s=idx(jj); yy=rtmed(s);
    h=errorbar(jj,yy,yy-rtq1(s),rtq3(s)-yy,'s','LineWidth',1.35,'MarkerSize',6);
    apply_label_style(h,labels7{jj});
end
grid on; xlim([0.35 numel(idx)+0.65]); set(gca,'XTick',1:numel(idx),'XTickLabel',labels7Tick,'XTickLabelRotation',28);
ylabel('Median complete runtime (s)'); if all(runtimeVals(isfinite(runtimeVals))>0); set(gca,'YScale','log'); end
title('(b) Executable runtime');
if ~isempty(maLocal) && isfinite(runtimeVals(maLocal)) && runtimeVals(maLocal)>0
    for jj=1:numel(runtimeVals)
        if isfinite(runtimeVals(jj)) && runtimeVals(jj)>0
            text(jj,runtimeVals(jj)*1.38,sprintf('%.1fx MA',runtimeVals(jj)/runtimeVals(maLocal)), ...
                'HorizontalAlignment','center','FontSize',11.5,'FontWeight','bold','Clipping','on');
        end
    end
end
yl=ylim; if all(yl>0); ylim([yl(1) yl(2)*1.55]); end

nexttile;
it=0:numel(D7.rad_gain_med)-1; plot_band(it,D7.rad_gain_med,D7.rad_gain_q1,D7.rad_gain_q3);
grid on;
stageBoundary=NaN;
if isfield(D7,'Stage1Accepts')
    z=D7.Stage1Accepts(isfinite(D7.Stage1Accepts));
    if ~isempty(z); stageBoundary=median(z); xline(stageBoundary,'--','HandleVisibility','off'); end
end
if isfinite(stageBoundary); text(0.03,0.14,sprintf('Stage-I end\nStage-II start'),'Units','normalized','FontSize',10.5,'FontWeight','bold','VerticalAlignment','bottom'); end
xlabel('Accepted update index'); ylabel('Radar-SNR gain (dB)');
title('(c) True-objective convergence');

nexttile; hold on;
gapBase=compact_labels(D7.labels(D7.baseline_idx));
gapLabels=[gapBase {'MA-RCG-MF (reference)'}];
gapTick=algorithm_tick_labels(gapLabels);
x=1:numel(gapLabels); y=[D7.GainMean_dB(:).' 0];
lo=[D7.GainLo_dB(:).' 0]; hi=[D7.GainHi_dB(:).' 0];
for jj=1:numel(x)
    h=errorbar(x(jj),y(jj),y(jj)-lo(jj),hi(jj)-y(jj),'d','LineWidth',1.35,'MarkerSize',6);
    apply_label_style(h,gapLabels{jj});
end
yline(0,'--k','HandleVisibility','off'); grid on;
set(gca,'XTick',x,'XTickLabel',gapTick,'XTickLabelRotation',28);
ylabel('MA-RCG-MF gain (dB)');
title('(d) Paired MA-RCG-MF gain (95% CI)');
text(0.02,0.95,'MA better','Units','normalized','FontSize',11,'FontWeight','bold','VerticalAlignment','top');
text(0.02,0.04,'Comparator better','Units','normalized','FontSize',11,'FontWeight','bold','VerticalAlignment','bottom');
text(0.98,0.52,'Equal','Units','normalized','FontSize',10.5,'HorizontalAlignment','right','VerticalAlignment','bottom');
sgtitle(tl,'Algorithm Performance, Runtime and Convergence','FontSize',18,'FontWeight','bold');
save_master(fig,output_dir,'F6_Algorithm_Validation'); close(fig);

fprintf('\nSix manuscript figures generated in:\n%s\n',output_dir);
fprintf('Experiment 12 is intentionally reported as complexity/significance tables, not a separate figure.\n');
end

function D=load_required(folder,name)
file=fullfile(folder,name); if ~exist(file,'file'); error('Required result is missing: %s',file); end
D=load(file);
end

function labels=nmse_tick_labels(v)
labels=cell(1,numel(v));
for ii=1:numel(v)
    if isinf(v(ii)) && v(ii)<0
        labels{ii}='Perfect';
    else
        labels{ii}=sprintf('%g',v(ii));
    end
end
end

function labels=numeric_tick_labels(v)
labels=arrayfun(@(x)sprintf('%g',x),v(:).','UniformOutput',false);
end

function labels=compact_labels(labels)
labels=cellstr(string(labels));
for i=1:numel(labels)
    % Replace exact long-form names before generic substrings.  In particular,
    % this prevents 'Proposed RET-RIS-RadarSNR' from becoming 'Proposed Proposed'.
    labels{i}=strrep(labels{i},'Proposed MA-RCG-MF','MA-RCG-MF');
    labels{i}=strrep(labels{i},'Accelerated Stage-I RCG-AO','Stage-I RCG-AO');
    labels{i}=strrep(labels{i},'Proposed RET-RIS-RadarSNR','MA-RCG-MF');
    labels{i}=strrep(labels{i},'Fixed RET, RIS-RadarSNR','Fixed RET/RIS');
    labels{i}=strrep(labels{i},'Fixed RET, No RIS','Fixed/no RIS');
    labels{i}=strrep(labels{i},'Sequential RET then RIS','Sequential');
    labels{i}=strrep(labels{i},'Fine AO block-coordinate phase search','Block Search');
    labels{i}=strrep(labels{i},'Fine AO block-coordinate search','Block Search');
    labels{i}=strrep(labels{i},'AO block-coordinate phase search','Block Search');
    labels{i}=strrep(labels{i},'AO block-coordinate search','Block Search');
    labels{i}=strrep(labels{i},'RET-RIS-RadarSNR','MA-RCG-MF');
    labels{i}=strrep(labels{i},'RET-RIS-MinPc','Joint RET-RIS MinPc');
    labels{i}=strrep(labels{i},'RET-NoRIS-MinPc','RET/no RIS');
    labels{i}=strrep(labels{i},'Fixed RET, RIS-MinPc','Fixed RET/RIS MinPc');
    labels{i}=strtrim(regexprep(labels{i},'\s+',' '));
end
end

function out=algorithm_tick_labels(labels)
% Short labels prevent x-tick text from crossing tile boundaries.
labels=cellstr(string(labels));
out=labels;
for ii=1:numel(labels)
    lab=lower(labels{ii});
    if contains(lab,'ma-rcg-mf') && contains(lab,'reference')
        out{ii}='MA ref.';
    elseif contains(lab,'ma-rcg-mf')
        out{ii}='MA-RCG-MF';
    elseif contains(lab,'block search')
        out{ii}='Block';
    elseif contains(lab,'stage-i')
        out{ii}='Stage-I';
    elseif contains(lab,'sequential')
        out{ii}='Sequential';
    elseif contains(lab,'fixed ret/ris') || contains(lab,'fixed/ris')
        out{ii}='Fixed RIS';
    elseif contains(lab,'fixed/no ris') || contains(lab,'fixed ret/no ris') || contains(lab,'no-ris')
        out{ii}='No RIS';
    else
        out{ii}=labels{ii};
    end
end
end

function apply_three_scheme_style(h,idx)
%APPLY_THREE_SCHEME_STYLE Consistent manuscript style for 3-scheme plots.
% 1: Fixed RET/no RIS    -> blue, solid, circle
% 2: RET/RIS baseline   -> orange, dashed, square
% 3: Proposed MA-RCG-MF -> yellow, dash-dot, diamond
colors=[0.0000 0.4470 0.7410; ...
        0.8500 0.3250 0.0980; ...
        0.9290 0.6940 0.1250];
lineStyles={'-','--','-.'};
markers={'o','s','d'};
idx=max(1,min(3,idx));
set(h,'Color',colors(idx,:), ...
    'LineStyle',lineStyles{idx}, ...
    'Marker',markers{idx}, ...
    'MarkerFaceColor','none');
end

function apply_label_style(h,label)
%APPLY_LABEL_STYLE Stable visual identity across all manuscript figures.
lab=lower(char(string(label)));
if contains(lab,'ma-rcg-mf')
    c=[0.9290 0.6940 0.1250]; ls='-.'; mk='d';
elseif contains(lab,'block search')
    c=[0.4940 0.1840 0.5560]; ls='-'; mk='^';
elseif contains(lab,'stage-i')
    c=[0.4660 0.6740 0.1880]; ls='--'; mk='v';
elseif contains(lab,'sequential')
    c=[0.3010 0.7450 0.9330]; ls=':'; mk='>';
elseif contains(lab,'joint ret-ris') || contains(lab,'joint minpc')
    c=[0.4940 0.1840 0.5560]; ls='-'; mk='p';
elseif contains(lab,'ret/no ris') || contains(lab,'ret only')
    c=[0.3010 0.7450 0.9330]; ls=':'; mk='^';
elseif contains(lab,'fixed ret/ris') || contains(lab,'fixed/ris') || contains(lab,'ris only')
    c=[0.8500 0.3250 0.0980]; ls='--'; mk='s';
elseif contains(lab,'fixed/no ris') || contains(lab,'fixed ret/no ris') || contains(lab,'no-ris')
    c=[0.0000 0.4470 0.7410]; ls='-'; mk='o';
else
    c=[0.25 0.25 0.25]; ls='-'; mk='o';
end
set(h,'Color',c,'LineStyle',ls,'Marker',mk,'MarkerFaceColor','none');
end

function apply_reference_style(h)
set(h,'Color',[0.35 0.35 0.35],'LineStyle','--','Marker','o', ...
    'MarkerFaceColor','none');
end

function set_feasibility_ylim(ax,Y)
% Zoom only when all observed rates are already high; the axis label remains
% explicit so the zoom is not mistaken for a zero-based bar chart.
v=Y(isfinite(Y));
if isempty(v); ylim(ax,[0 105]); return; end
mn=min(v);
if mn>=60
    lo=max(0,10*floor((mn-10)/10));
    ylim(ax,[lo 102]);
else
    ylim(ax,[0 102]);
end
end

function annotate_heatmap_values(ax,A,climmax)
% Print numerical values so the heat map remains interpretable in grayscale.
for rr=1:size(A,1)
    for cc=1:size(A,2)
        v=A(rr,cc); if ~isfinite(v); continue; end
        if climmax>0 && v/climmax>0.58
            tc=[0 0 0];
        else
            tc=[1 1 1];
        end
        text(ax,cc,rr,sprintf('%.2f',v),'HorizontalAlignment','center', ...
            'VerticalAlignment','middle','Color',tc,'FontSize',13.5,'FontWeight','bold');
    end
end
end

function bar_ci(mu,lo,hi,labels)
mu=mu(:).'; lo=lo(:).'; hi=hi(:).'; x=1:numel(mu);
b=bar(x,mu,'FaceColor','flat'); hold on;
for ii=1:numel(labels)
    b.CData(ii,:)=scheme_color(labels{ii});
end
errorbar(x,mu,mu-lo,hi-mu,'.k','LineWidth',1.1);
set(gca,'XTick',x,'XTickLabel',labels,'XTickLabelRotation',25);
end

function c=scheme_color(label)
lab=lower(char(string(label)));
if contains(lab,'ma-rcg-mf')
    c=[0.9290 0.6940 0.1250];
elseif contains(lab,'joint')
    c=[0.4940 0.1840 0.5560];
elseif contains(lab,'ret only') || contains(lab,'ret/no ris')
    c=[0.3010 0.7450 0.9330];
elseif contains(lab,'ris only') || contains(lab,'fixed/ris') || contains(lab,'fixed ret/ris')
    c=[0.8500 0.3250 0.0980];
else
    c=[0.0000 0.4470 0.7410];
end
end

function [med,q1,q3]=runtime_iqr(X)
med=nan(1,size(X,2)); q1=med; q3=med;
for j=1:size(X,2)
    v=X(:,j); v=v(isfinite(v));
    if ~isempty(v); med(j)=median(v); q1(j)=empirical_quantile(v,0.25); q3(j)=empirical_quantile(v,0.75); end
end
end

function plot_band(x,med,q1,q3)
x=x(:).'; med=med(:).'; q1=q1(:).'; q3=q3(:).'; ok=isfinite(x)&isfinite(med)&isfinite(q1)&isfinite(q3);
if ~any(ok); return; end
xx=x(ok); mm=med(ok); ll=q1(ok); uu=q3(ok);
fill([xx fliplr(xx)],[ll fliplr(uu)],[0.85 0.85 0.85], ...
    'EdgeColor','none','FaceAlpha',0.45); hold on;
plot(xx,mm,'-o','LineWidth',1.5,'MarkerSize',6);
end

function save_master(fig,folder,name)
set(findall(fig,'-property','FontName'),'FontName','Times New Roman');
% Collision-safe typography: large enough for two-column placement, but
% annotations are NOT globally enlarged because that can make labels collide.
axs0=findall(fig,'Type','axes');
for kk0=1:numel(axs0)
    try
        axs0(kk0).FontSize=14.5;
        axs0(kk0).LineWidth=1.15;
        axs0(kk0).TickLength=[0.014 0.014];
        axs0(kk0).LooseInset=max(axs0(kk0).TightInset,[0.02 0.02 0.02 0.02]);
    catch
    end
    try
        axs0(kk0).XLabel.FontSize=16.0;
        axs0(kk0).YLabel.FontSize=16.0;
        axs0(kk0).XLabel.FontWeight='normal';
        axs0(kk0).YLabel.FontWeight='normal';
    catch
    end
    try
        axs0(kk0).Title.FontSize=15.5;
        axs0(kk0).Title.FontWeight='normal';
    catch
    end
end
try; set(findall(fig,'Type','legend'),'FontSize',12.5,'Box','off'); catch; end
try; set(findall(fig,'Type','ColorBar'),'FontSize',13.0); catch; end
% Preserve explicitly chosen annotation sizes instead of forcing every text
% object to >=14 pt; that was the main cause of several collisions.
ln0=findall(fig,'Type','line');
for kk0=1:numel(ln0)
    try; if ln0(kk0).LineWidth < 1.8; ln0(kk0).LineWidth=1.8; end; catch; end
    try; if ~strcmp(ln0(kk0).Marker,'none') && ln0(kk0).MarkerSize < 7; ln0(kk0).MarkerSize=7; end; catch; end
end
axs=findall(fig,'Type','axes');
for kk=1:numel(axs)
    try; axs(kk).Toolbar.Visible='off'; catch; end
end
drawnow;
png=fullfile(folder,[name '.png']); pdf=fullfile(folder,[name '.pdf']);
figfile=fullfile(folder,[name '.fig']);
try; exportgraphics(fig,png,'Resolution',600); catch; print(fig,png,'-dpng','-r600'); end
try; exportgraphics(fig,pdf,'ContentType','vector'); catch; print(fig,pdf,'-dpdf','-painters'); end
try; savefig(fig,figfile); catch ME; warning('Could not save %s: %s',figfile,ME.message); end
end
