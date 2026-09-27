% EXPORT_ALL_PAPER_STATISTICS
% Export the exact PAPER-mode statistics produced by the statistically
% final paper experiments experiments. Run after Experiments 1--12.
clear; clc;
root=fileparts(mfilename('fullpath')); if isempty(root); root=pwd; end
addpath(root); addpath(fullfile(root,'function'));
results_dir=fullfile(root,'results');
out_dir=fullfile(results_dir,'paper_tables');
if ~exist(out_dir,'dir'); mkdir(out_dir); end

scheme4={'Fixed RET, No RIS','Fixed RET, RIS-MinPc', ...
    'RET-MinPc, No RIS','RET-RIS-MinPc'};
scheme5={'Fixed RET, No RIS','Fixed RET, RIS-MinPc', ...
    'RET-MinPc, No RIS','RET-RIS-MinPc','RET-RIS-RadarSNR'};
scheme3={'Fixed RET, No RIS','RET-RIS-MinPc','RET-RIS-RadarSNR'};
missing={};

%% Experiment 1: already contains the exact common-feasible and paired tables.
f=fullfile(results_dir,'exp1_revised_paper.mat');
if exist(f,'file')
    D=load(f);
    if isfield(D,'Tpower'); writetable(D.Tpower,fullfile(out_dir,'Exp1_Communication_Power_Common_Feasible.csv')); end
    if isfield(D,'Tgain'); writetable(D.Tgain,fullfile(out_dir,'Exp1_Paired_Marginal_Gains.csv')); end
else
    missing{end+1}=f; %#ok<AGROW>
end

%% Experiment 2
f=fullfile(results_dir,'exp2_revised_paper.mat');
if exist(f,'file')
    D=load(f);
    writetable(curve_common_table(D.range_P_dBm,'Ptotal_dBm',D.Obj,D.Feasible, ...
        D.CommonMain,scheme5,'RadarSNR','dB'),fullfile(out_dir,'Exp2_Radar_SNR_Common_Feasible.csv'));
    writetable(curve_common_table(D.range_P_dBm,'Ptotal_dBm',D.Pcomm,D.Feasible, ...
        D.CommonMain,scheme5,'Pcomm','dBm'),fullfile(out_dir,'Exp2_Communication_Power_Common_Feasible.csv'));
    writetable(curve_common_table(D.range_P_dBm,'Ptotal_dBm',D.Pradar,D.Feasible, ...
        D.CommonMain,scheme5,'Pradar','W'),fullfile(out_dir,'Exp2_Residual_Radar_Power_Common_Feasible.csv'));

    allocation_labels={'Dynamic fixed RET/no RIS','Dynamic RET-RIS-RadarSNR', ...
        'Fixed radar fraction 0.5','Fixed radar fraction 0.7'};
    A=cat(3,D.Obj(:,:,1),D.Obj(:,:,5));
    if isfield(D,'ObjAllocation_dB')
        % Reconstruct the two displayed fixed-split samples from the raw arrays.
        [~,a50]=min(abs(D.cfg.alpha_fixed_list-0.5));
        [~,a70]=min(abs(D.cfg.alpha_fixed_list-0.7));
        A=cat(3,A,D.Obj_fixed_alpha(:,:,a50),D.Obj_fixed_alpha(:,:,a70));
        F=cat(3,D.Feasible(:,:,1),D.Feasible(:,:,5), ...
            D.Fea_fixed_alpha(:,:,a50),D.Fea_fixed_alpha(:,:,a70));
        writetable(curve_common_table(D.range_P_dBm,'Ptotal_dBm',A,F, ...
            D.CommonAllocation,allocation_labels,'RadarSNR','dB'), ...
            fullfile(out_dir,'Exp2_Full_Budget_Fixed_Split_Comparison.csv'));
    end

    if isfield(D,'Pmin_samples_W')
        writetable(inverse_power_table(D.radar_targets_dB,D.Pmin_samples_W, ...
            D.Pmin_common_mask,scheme5), ...
            fullfile(out_dir,'Exp2_Equal_Sensing_Minimum_Power.csv'));
    end
else
    missing{end+1}=f; %#ok<AGROW>
end

%% Experiment 3
f=fullfile(results_dir,'exp3_qos_revised_paper.mat');
if exist(f,'file')
    D=load(f);
    writetable(curve_common_table(D.gamma_dB_range,'SINR_dB',D.Obj,D.Feasible, ...
        D.CommonFeasible,scheme5,'RadarSNR','dB'),fullfile(out_dir,'Exp3_Radar_SNR_Common_Feasible.csv'));
    writetable(curve_common_table(D.gamma_dB_range,'SINR_dB',D.Pcomm,D.Feasible, ...
        D.CommonFeasible,scheme5,'Pcomm','dBm'),fullfile(out_dir,'Exp3_Communication_Power_Common_Feasible.csv'));
    writetable(curve_common_table(D.gamma_dB_range,'SINR_dB',D.Pradar,D.Feasible, ...
        D.CommonFeasible,scheme5,'Pradar','W'),fullfile(out_dir,'Exp3_Residual_Radar_Power_Common_Feasible.csv'));
    if isfield(D,'RadarGain_dB')
        T=array2table([D.gamma_dB_range(:),D.RadarGain_dB(:,2:5), ...
            D.RadarGain_lo(:,2:5),D.RadarGain_hi(:,2:5),D.RadarGain_n(:,2:5)], ...
            'VariableNames',{'SINR_dB','Gain2_dB','Gain3_dB','Gain4_dB','Gain5_dB', ...
            'Gain2_CI95_Low','Gain3_CI95_Low','Gain4_CI95_Low','Gain5_CI95_Low', ...
            'Gain2_CI95_High','Gain3_CI95_High','Gain4_CI95_High','Gain5_CI95_High', ...
            'Gain2_N','Gain3_N','Gain4_N','Gain5_N'});
        writetable(T,fullfile(out_dir,'Exp3_Paired_Radar_SNR_Gains.csv'));
    end
else
    missing{end+1}=f; %#ok<AGROW>
end

%% Experiment 4
f=fullfile(results_dir,'exp4_ablation_revised_paper.mat');
if exist(f,'file')
    D=load(f);
    if isfield(D,'Tsummary'); writetable(D.Tsummary,fullfile(out_dir,'Exp4_Dual_Mechanism_Ablation.csv')); end
    if isfield(D,'Tsame'); writetable(D.Tsame,fullfile(out_dir,'Exp4_Same_Objective_Hardware_Ablation.csv')); end
else
    missing{end+1}=f; %#ok<AGROW>
end

%% Experiment 5
f=fullfile(results_dir,'exp5_scalability_revised_paper.mat');
if exist(f,'file')
    D=load(f);
    writetable(curve_common_table(D.N_range,'RIS_elements',D.ObjN,D.FeaN, ...
        D.CommonN,scheme3,'RadarSNR','dB'),fullfile(out_dir,'Exp5_Radar_SNR_vs_N_Common_Feasible.csv'));
    writetable(curve_common_table(D.K_range,'Users',D.ObjK,D.FeaK, ...
        D.CommonK,scheme3,'RadarSNR','dB'),fullfile(out_dir,'Exp5_Radar_SNR_vs_K_Common_Feasible.csv'));
    writetable(curve_common_table(D.N_range,'RIS_elements',D.PcommN,D.FeaN, ...
        D.CommonN,scheme3,'Pcomm','dBm'),fullfile(out_dir,'Exp5_Pcomm_vs_N_Common_Feasible.csv'));
    writetable(curve_common_table(D.K_range,'Users',D.PcommK,D.FeaK, ...
        D.CommonK,scheme3,'Pcomm','dBm'),fullfile(out_dir,'Exp5_Pcomm_vs_K_Common_Feasible.csv'));
    if isfield(D,'ObjNGain_dB')
        T=array2table([D.N_range(:),D.ObjNGain_dB,D.ObjNGain_lo,D.ObjNGain_hi,D.ObjNGain_n], ...
            'VariableNames',{'RIS_elements','Gain1_dB','Gain2_dB','Gain3_dB', ...
            'Gain1_CI95_Low','Gain2_CI95_Low','Gain3_CI95_Low', ...
            'Gain1_CI95_High','Gain2_CI95_High','Gain3_CI95_High', ...
            'Gain1_N','Gain2_N','Gain3_N'});
        writetable(T,fullfile(out_dir,'Exp5_Paired_RIS_Size_Gain.csv'));
    end
    if isfield(D,'TimeN_med') && any(isfinite(D.TimeN_med(:)))
        writetable(timing_table(D.N_range,'RIS_elements',D.TimeN_med,D.TimeN_q1,D.TimeN_q3,scheme3), ...
            fullfile(out_dir,'Exp5_Runtime_vs_N.csv'));
        writetable(timing_table(D.K_range,'Users',D.TimeK_med,D.TimeK_q1,D.TimeK_q3,scheme3), ...
            fullfile(out_dir,'Exp5_Runtime_vs_K.csv'));
    end
else
    missing{end+1}=f; %#ok<AGROW>
end

%% Experiment 6
f=fullfile(results_dir,'exp6_multitimescale_mobility_paper.mat');
if exist(f,'file')
    D=load(f); rows=cell(0,20);
    for vi=1:numel(D.speed_kmh)
        for ti=1:numel(D.update_ms)
            x=D.Loss_dB(:,vi,ti); x=x(isfinite(x));
            [mu,lo,hi,n]=scalar_ci(x);
            cp=D.CommChange_dB(:,vi,ti); cp=cp(isfinite(cp));
            [cp_mu,cp_lo,cp_hi,cp_n]=scalar_ci(cp);
            fallback_count=sum(D.FreshDominanceFallback(:,vi,ti));
            rows(end+1,:)={D.speed_kmh(vi),D.update_ms(ti), ...
                (D.speed_kmh(vi)/3.6)*(D.update_ms(ti)/1000), ...
                mu,lo,hi,D.LossP90(vi,ti),D.LossP90Lo(vi,ti),D.LossP90Hi(vi,ti), ...
                D.LossP90Count(vi,ti),n,cp_mu,cp_lo,cp_hi,cp_n, ...
                D.StaleRate(vi,ti),D.FreshRate(vi,ti),D.StaleOutage(vi,ti), ...
                D.StressOutage(vi,ti),fallback_count}; %#ok<AGROW>
        end
    end
    T=cell2table(rows,'VariableNames',{'Speed_kmh','RIS_update_ms','Displacement_m', ...
        'Loss_mean_dB','Loss_CI95_low_dB','Loss_CI95_high_dB','Loss_P90_dB', ...
        'Loss_P90_bootstrap_CI95_low_dB','Loss_P90_bootstrap_CI95_high_dB', ...
        'Loss_P90_samples','Loss_mean_samples','CommChange_mean_dB', ...
        'CommChange_CI95_low_dB','CommChange_CI95_high_dB','CommChange_samples', ...
        'Stale_feasibility_percent','Fresh_feasibility_percent', ...
        'Conditional_outage_percent','Stress_outage_percent', ...
        'Fresh_dominance_fallback_count'});
    writetable(T,fullfile(out_dir,'Exp6_Mobility_Aging.csv'));
    if isfield(D,'Ttrend') && istable(D.Ttrend)
        writetable(D.Ttrend,fullfile(out_dir,'Exp6_Paired_Trend_Audit.csv'));
    end
    srcAudit=fullfile(results_dir,'exp6_pairing_audit_paper.txt');
    if exist(srcAudit,'file'); copyfile(srcAudit,fullfile(out_dir,'Exp6_Pairing_Audit.txt')); end
else
    missing{end+1}=f; %#ok<AGROW>
end

%% Experiment 7
f=fullfile(results_dir,'exp7_algorithm_convergence_paper.mat');
if exist(f,'file')
    D=load(f);
    if isfield(D,'Tsummary'); writetable(D.Tsummary,fullfile(out_dir,'Exp7_Algorithm_Comparison.csv')); end
    if isfield(D,'Tpaired'); writetable(D.Tpaired,fullfile(out_dir,'Exp7_Paired_Objective_Gaps_Common_Main.csv')); end
    if isfield(D,'Trefine'); writetable(D.Trefine,fullfile(out_dir,'Exp7_Enhanced_Refinement_Per_MC.csv')); end
    if isfield(D,'Truntime'); writetable(D.Truntime,fullfile(out_dir,'Exp7_Runtime_Breakdown_Per_MC.csv')); end
    if isfield(D,'Tv2'); writetable(D.Tv2,fullfile(out_dir,'Exp7_V2_StageII_Metrics_Per_MC.csv')); end
    audit_src=fullfile(results_dir,'exp7_v2_postrun_audit_paper.csv');
    if exist(audit_src,'file'); copyfile(audit_src,fullfile(out_dir,'Exp7_V2_Postrun_Audit.csv')); end
else
    missing{end+1}=f; %#ok<AGROW>
end


%% Experiments 8--12: copy their already-structured PAPER CSV tables.
extra_files={ ...
    'exp8_ris_deployment_paper.csv','Exp8_RIS_Deployment.csv'; ...
    'exp9_csi_uncertainty_paper.csv','Exp9_CSI_Sensitivity.csv'; ...
    'exp10_ris_quantization_paper.csv','Exp10_RIS_Post_Quantization.csv'; ...
    'exp11_ret_resolution_paper.csv','Exp11_RET_Resolution.csv'; ...
    'exp12_complexity_paper.csv','Exp12_Complexity.csv'; ...
    'exp12_performance_paper.csv','Exp12_Performance.csv'; ...
    'exp12_significance_paper.csv','Exp12_Significance.csv'; ...
    'exp12_solver_validation_paper.csv','Exp12_Solver_Validation.csv'; ...
    'exp12_solver_validation_summary_paper.csv','Exp12_Solver_Validation_Summary.csv'};
for ii=1:size(extra_files,1)
    src=fullfile(results_dir,extra_files{ii,1});
    if exist(src,'file')
        copyfile(src,fullfile(out_dir,extra_files{ii,2}));
    else
        missing{end+1}=src; %#ok<AGROW>
    end
end

fprintf('\nFinal PAPER tables exported to:\n%s\n',out_dir);
if ~isempty(missing)
    fprintf('\nThe following PAPER-mode MAT files were not found:\n');
    for ii=1:numel(missing); fprintf('  %s\n',missing{ii}); end
    fprintf('Run the corresponding experiment(s) in PAPER mode, then rerun this exporter.\n');
end

%% Local helpers
function T=curve_common_table(axis_values,axis_name,X,F,common_mask,labels,metric,unit)
if ndims(X)==2; X=reshape(X,size(X,1),size(X,2),1); end
if ndims(F)==2; F=reshape(F,size(F,1),size(F,2),1); end
[I,P,S]=size(X);
if ~isequal(size(common_mask),[I P]); error('Common mask size mismatch.'); end
rows=cell(0,11);
for pp=1:P
    for ss=1:S
        ok=common_mask(:,pp)&F(:,pp,ss)&isfinite(X(:,pp,ss));
        v=X(ok,pp,ss);
        [mu,lo,hi,n]=metric_ci(v,unit);
        rows(end+1,:)={axis_values(pp),labels{ss},metric,unit,mu,lo,hi,n, ...
            100*sum(common_mask(:,pp))/I,100*mean(F(:,pp,ss)),I}; %#ok<AGROW>
    end
end
T=cell2table(rows,'VariableNames',{matlab.lang.makeValidName(axis_name), ...
    'Scheme','Metric','Unit','Mean','CI95_Low','CI95_High','Common_samples', ...
    'Common_feasible_percent','Scheme_feasible_percent','Total_samples'});
end

function T=inverse_power_table(targets,Pmin,common_mask,labels)
[I,Q,S]=size(Pmin); rows=cell(0,10);
for qq=1:Q
    for ss=1:S
        ok=common_mask(:,qq)&isfinite(Pmin(:,qq,ss));
        [mu,lo,hi,n]=metric_ci(Pmin(ok,qq,ss),'dBm');
        rows(end+1,:)={targets(qq),labels{ss},mu,lo,hi,n, ...
            100*sum(common_mask(:,qq))/I,100*mean(isfinite(Pmin(:,qq,ss))),I,'dBm'}; %#ok<AGROW>
    end
end
T=cell2table(rows,'VariableNames',{'RadarSNR_target_dB','Scheme','MinimumPower_mean_dBm', ...
    'CI95_Low_dBm','CI95_High_dBm','Common_samples','Common_reachable_percent', ...
    'Scheme_reachable_percent','Total_samples','Unit'});
end

function T=timing_table(axis_values,axis_name,med,q1,q3,labels)
rows=cell(0,6);
for pp=1:numel(axis_values)
    for ss=1:numel(labels)
        rows(end+1,:)={axis_values(pp),labels{ss},med(pp,ss),q1(pp,ss),q3(pp,ss),'s'}; %#ok<AGROW>
    end
end
T=cell2table(rows,'VariableNames',{matlab.lang.makeValidName(axis_name), ...
    'Scheme','Median','Q1','Q3','Unit'});
end

function [mu,lo,hi,n]=metric_ci(x,unit)
x=x(isfinite(x)); n=numel(x);
if n==0; mu=NaN; lo=NaN; hi=NaN; return; end
m=mean(x); if n>1; half=student_t_critical(0.95,n-1)*std(x,0)/sqrt(n); else; half=0; end
switch unit
    case 'dB'
        mu=10*log10(max(m,realmin));
        lo=10*log10(max(m-half,realmin));
        hi=10*log10(max(m+half,realmin));
    case 'dBm'
        mu=10*log10(max(m,realmin)/1e-3);
        lo=10*log10(max(m-half,realmin)/1e-3);
        hi=10*log10(max(m+half,realmin)/1e-3);
    otherwise
        mu=m; lo=m-half; hi=m+half;
end
end

function [mu,lo,hi,n]=scalar_ci(x)
x=x(isfinite(x)); n=numel(x);
if n==0; mu=NaN; lo=NaN; hi=NaN; return; end
mu=mean(x); if n>1; d=student_t_critical(0.95,n-1)*std(x,0)/sqrt(n); else; d=0; end
lo=mu-d; hi=mu+d;
end
