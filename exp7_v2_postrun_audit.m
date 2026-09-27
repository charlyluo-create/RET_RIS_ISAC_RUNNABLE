% =========================================================================
% EXP7_V2_POSTRUN_AUDIT
% Post-run audit for MA-RCG-MF v2 Stage-II.
%
% Run AFTER exp7_algorithm_comparison_convergence_final.m. This script does
% not rerun any optimizer. It checks whether the diagnosis-driven V2 closes
% the Stage-II objective gap while preserving a meaningful exact-evaluation
% and runtime advantage over the fine Block Search reference.
% =========================================================================
clear; close all; clc;

RUN_MODE=getenv('RET_RIS_RUN_MODE');
if isempty(RUN_MODE); RUN_MODE='paper'; end
RUN_MODE=lower(char(RUN_MODE));
root=fileparts(mfilename('fullpath')); if isempty(root); root=pwd; end
addpath(root); addpath(fullfile(root,'function'));
cfg=unified_ret_ris_isac_config(RUN_MODE);
results_dir=fullfile(root,'results');
source_file=fullfile(results_dir,sprintf('exp7_algorithm_convergence_%s.mat',RUN_MODE));
if exist(source_file,'file')~=2
    error('Missing %s. Run Experiment 7 first.',source_file);
end
D=load(source_file);
if ~isfield(D,'cfg') || ~isfield(D.cfg,'paper_release_version') || ...
        ~strcmp(D.cfg.paper_release_version,cfg.paper_release_version)
    error('Exp.7 MAT belongs to a different release. Re-run V2 Exp.7.');
end
required={'Obj','Feasible','RefineGain_dB','RefineExactEvals','BCSExactEvals', ...
    'Runtime','Stage2Runtime','BCSRefineRuntime','V2Metrics'};
for ii=1:numel(required)
    if ~isfield(D,required{ii}); error('Missing Exp.7 field: %s',required{ii}); end
end

stage1=7; block=5; ma=8;
common=D.Feasible(:,stage1)&D.Feasible(:,block)&D.Feasible(:,ma)& ...
    isfinite(D.Obj(:,stage1))&isfinite(D.Obj(:,block))&isfinite(D.Obj(:,ma))& ...
    D.Obj(:,stage1)>0&D.Obj(:,block)>0&D.Obj(:,ma)>0;
if nnz(common)<1; error('No common-feasible Stage-I/Block/MA samples.'); end
stage1_db=10*log10(D.Obj(common,stage1));
ma_db=10*log10(D.Obj(common,ma));
block_db=10*log10(D.Obj(common,block));
ma_gain=ma_db-stage1_db;
block_gain=block_db-stage1_db;
ma_minus_block=ma_db-block_db;

ma_exact=D.RefineExactEvals(common);
bcs_exact=D.BCSExactEvals(common);
ma_rt=D.Runtime(common,ma);
bcs_rt=D.Runtime(common,block);
stage2_rt=D.Stage2Runtime(common);
bcs_ref_rt=D.BCSRefineRuntime(common);
V=D.V2Metrics(common,:);

stats=@(x)[mean_finite(x),median_finite(x),qfinite(x,0.25),qfinite(x,0.75)];
A=stats(ma_gain); B=stats(block_gain); G=stats(ma_minus_block);
ER=stats(bcs_exact./max(ma_exact,eps)); RR=stats(bcs_rt./max(ma_rt,eps));
RRef=stats(bcs_ref_rt./max(stage2_rt,eps));

fprintf('\n===============================================================\n');
fprintf(' MA-RCG-MF v2 POST-RUN AUDIT\n');
fprintf(' Release: %s\n',cfg.paper_release_version);
fprintf(' Common samples: %d/%d\n',nnz(common),size(D.Obj,1));
fprintf('===============================================================\n');
fprintf('Stage-II gain, MA v2: mean %.4f dB; median %.4f dB\n',A(1),A(2));
fprintf('Stage-II gain, Block: mean %.4f dB; median %.4f dB\n',B(1),B(2));
fprintf('MA v2 minus Block:   mean %.4f dB; median %.4f dB\n',G(1),G(2));
fprintf('Block/MA exact-eval ratio: median %.2fx\n',ER(2));
fprintf('Block/MA complete-runtime ratio: median %.2fx\n',RR(2));
fprintf('Block/MA refinement-runtime ratio: median %.2fx\n',RRef(2));
fprintf('V2 sweeps completed: median %.1f\n',median_finite(V(:,1)));
fprintf('V2 block visits: median %.1f\n',median_finite(V(:,2)));
fprintf('V2 budget expansions: median %.1f\n',median_finite(V(:,3)));
fprintf('V2 RIS exact evaluations: median %.1f\n',median_finite(V(:,5)));
fprintf('V2 final mean trust radius: median %.2f deg\n',median_finite(V(:,6)));
fprintf('V2 accepted block updates: median %.1f\n',median_finite(V(:,8)));

report_file=fullfile(results_dir,sprintf('exp7_v2_postrun_audit_%s.txt',RUN_MODE));
fid=fopen(report_file,'w'); if fid<0; error('Cannot write %s.',report_file); end
cleaner=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'MA-RCG-MF V2 POST-RUN AUDIT\n');
fprintf(fid,'Release: %s\n',cfg.paper_release_version);
fprintf(fid,'Common samples: %d/%d\n\n',nnz(common),size(D.Obj,1));
fprintf(fid,'MA v2 Stage-II gain mean/median: %.6f / %.6f dB\n',A(1),A(2));
fprintf(fid,'Block Stage-II gain mean/median: %.6f / %.6f dB\n',B(1),B(2));
fprintf(fid,'MA v2 minus Block mean/median: %.6f / %.6f dB\n',G(1),G(2));
fprintf(fid,'Block/MA exact-eval ratio median: %.6f\n',ER(2));
fprintf(fid,'Block/MA complete-runtime ratio median: %.6f\n',RR(2));
fprintf(fid,'Block/MA refinement-runtime ratio median: %.6f\n',RRef(2));
fprintf(fid,'V2 sweeps median: %.6f\n',median_finite(V(:,1)));
fprintf(fid,'V2 block visits median: %.6f\n',median_finite(V(:,2)));
fprintf(fid,'V2 budget expansions median: %.6f\n',median_finite(V(:,3)));
fprintf(fid,'V2 RET exact evals median: %.6f\n',median_finite(V(:,4)));
fprintf(fid,'V2 RIS exact evals median: %.6f\n',median_finite(V(:,5)));
fprintf(fid,'V2 final mean/max trust median: %.6f / %.6f deg\n', ...
    median_finite(V(:,6)),median_finite(V(:,7)));
fprintf(fid,'V2 accepted block updates median: %.6f\n\n',median_finite(V(:,8)));
fprintf(fid,'Interpretation guide:\n');
fprintf(fid,'- |MA-Block| <= 0.25 dB with a clear runtime/evaluation advantage: strong target.\n');
fprintf(fid,'- MA-Block still below -0.4 dB: inspect per-block proxy ranking or increase sequential sweep capacity.\n');
fprintf(fid,'- Block/MA runtime ratio near 1 despite fewer exact solves: proxy/RCG overhead dominates.\n');
fprintf(fid,'- This audit never edits or smooths numerical results.\n');

T=table(find(common),ma_gain,block_gain,ma_minus_block,ma_exact,bcs_exact,ma_rt,bcs_rt, ...
    V(:,1),V(:,2),V(:,3),V(:,4),V(:,5),V(:,6),V(:,7),V(:,8), ...
    'VariableNames',{'MC_index','MA_StageII_gain_dB','Block_gain_dB','MA_minus_Block_dB', ...
    'MA_refine_exact_evals','Block_refine_exact_evals','MA_complete_runtime_s', ...
    'Block_complete_runtime_s','V2_sweeps','V2_block_visits','V2_budget_expansions', ...
    'V2_RET_exact','V2_RIS_exact','V2_mean_final_trust_deg','V2_max_final_trust_deg', ...
    'V2_accepted_block_updates'});
audit_csv=fullfile(results_dir,sprintf('exp7_v2_postrun_audit_%s.csv',RUN_MODE));
writetable(T,audit_csv);
% Close the report before packaging it for an easy chat handoff.
clear cleaner
send_zip=fullfile(results_dir,sprintf('EXP7_V2_SEND_TO_CHAT_%s.zip',upper(RUN_MODE)));
if exist(send_zip,'file'); delete(send_zip); end
handoff={report_file,audit_csv, ...
    fullfile(results_dir,sprintf('exp7_algorithm_summary_%s.csv',RUN_MODE)), ...
    fullfile(results_dir,sprintf('exp7_paired_objective_gaps_%s.csv',RUN_MODE)), ...
    fullfile(results_dir,sprintf('exp7_ma_rcg_mf_refinement_%s.csv',RUN_MODE)), ...
    fullfile(results_dir,sprintf('exp7_runtime_breakdown_%s.csv',RUN_MODE)), ...
    fullfile(results_dir,sprintf('exp7_v2_stage2_metrics_%s.csv',RUN_MODE))};
handoff=handoff(cellfun(@(f)exist(f,'file')==2,handoff));
if ~isempty(handoff)
    handoff_names=cell(size(handoff));
    for ii=1:numel(handoff)
        [~,nn,ee]=fileparts(handoff{ii}); handoff_names{ii}=[nn ee];
    end
    zip(send_zip,handoff_names,results_dir);
end
fprintf('\nAudit report: %s\n',report_file);
fprintf('Chat handoff ZIP: %s\n',send_zip);

function m=mean_finite(x)
x=x(isfinite(x)); if isempty(x); m=NaN; else; m=mean(x); end
end
function m=median_finite(x)
x=x(isfinite(x)); if isempty(x); m=NaN; else; m=median(x); end
end
function q=qfinite(x,p)
x=x(isfinite(x)); if isempty(x); q=NaN; else; q=empirical_quantile(x,p); end
end
