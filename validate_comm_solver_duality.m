function report = validate_comm_solver_duality()
%VALIDATE_COMM_SOLVER_DUALITY Cross-check FPI against the CVX/SOCP backend.
%
% The test intentionally spans several M/K/QoS operating points. For every
% case it verifies (i) feasibility on the original, unnormalized channel,
% (ii) minimum-power agreement with CVX, and (iii) the internal duality
% certificate. Results are saved so Experiment 12 can report solver accuracy
% and micro-benchmark speedup together with the complexity table.

root=fileparts(mfilename('fullpath')); if isempty(root); root=pwd; end
addpath(root); addpath(fullfile(root,'function'));
if exist('cvx_begin','file')~=2
    error('CVX is required for the duality-vs-CVX validation.');
end
mode=getenv('RET_RIS_RUN_MODE'); if isempty(mode); mode='fast'; end; mode=lower(mode);
cfg=unified_ret_ris_isac_config(mode);
results_dir=fullfile(root,'results'); if ~exist(results_dir,'dir'); mkdir(results_dir); end

% [M K gamma_dB seedOffset]
cases=[ ...
     8 3  5 1; ...
     8 4  8 2; ...
    12 4  8 3; ...
    12 4 12 4; ...
    16 4 10 5; ...
    16 4 12 6; ...
    16 4 16 7; ...
    16 6 10 8];
sigma=3.2e-6;
C=size(cases,1);
rel_gap=nan(C,1); t_fast=rel_gap; t_cvx=rel_gap; min_ratio=rel_gap;
fpi_iter=rel_gap; pd_gap=rel_gap; hybrid_gap=rel_gap; hybrid_backend=strings(C,1);

comm_solver_telemetry('reset');
for cc=1:C
    M=cases(cc,1); K=cases(cc,2); gdB=cases(cc,3); gamma=10^(gdB/10);
    rng(88000+cases(cc,4),'twister');
    H=(randn(M,K)+1i*randn(M,K))/sqrt(2);
    % A controlled common component prevents the test from being unrealistically
    % orthogonal while retaining full rank in the tested M>=K regime.
    H=H+0.12*repmat(H(:,1),1,K);
    H=2e-4*H/max(norm(H,'fro'),eps)*sqrt(M*K);

    tf=tic; [Wf,pf,ff,sf,mf]=solve_comm_min_power(H,gamma,sigma,struct('mode','duality')); t_fast(cc)=toc(tf);
    tc=tic; [~,pc,fc,sc]=solve_comm_min_power(H,gamma,sigma,struct('mode','cvx')); t_cvx(cc)=toc(tc);
    if ~ff; error('Duality solver failed validation case %d: %s',cc,sf); end
    if ~fc; error('CVX solver failed validation case %d: %s',cc,sc); end
    rel_gap(cc)=abs(pf-pc)/max(pc,1e-12);
    sinr=compute_sinr(H,Wf,sigma); min_ratio(cc)=min(sinr/gamma);
    fpi_iter(cc)=mf.fpi_iterations; pd_gap(cc)=mf.primal_dual_gap;
    if rel_gap(cc)>cfg.comm_solver_validation_power_gap_tol
        error('Duality/CVX power mismatch %.3g exceeds tolerance in case %d.',rel_gap(cc),cc);
    end
    if min_ratio(cc)<1-cfg.comm_solver_validation_sinr_tol
        error('Duality solver SINR validation failed in case %d.',cc);
    end
    if ~isfinite(pd_gap(cc)) || pd_gap(cc)>2e-4
        error('Duality primal-dual certificate failed in case %d.',cc);
    end

    % Verify the production hybrid path produces the same minimum power.
    [~,ph,fh,sh,mh]=solve_comm_min_power(H,gamma,sigma,struct('mode','hybrid'));
    if ~fh; error('Hybrid solver failed validation case %d: %s',cc,sh); end
    hybrid_gap(cc)=abs(ph-pc)/max(pc,1e-12);
    hybrid_backend(cc)=string(mh.backend_used);
    if hybrid_gap(cc)>cfg.comm_solver_validation_power_gap_tol
        error('Hybrid/CVX power mismatch %.3g exceeds tolerance in case %d.',hybrid_gap(cc),cc);
    end
end

speedup=t_cvx./max(t_fast,eps);
report=table((1:C).',cases(:,1),cases(:,2),cases(:,3),rel_gap,hybrid_gap, ...
    min_ratio,pd_gap,fpi_iter,t_fast,t_cvx,speedup,hybrid_backend, ...
    'VariableNames',{'Case','M','K','Gamma_dB','RelativePowerGap', ...
    'HybridRelativePowerGap','MinSINRRatio','PrimalDualGap','FPIIterations', ...
    'DualityTime_s','CVXTime_s','Speedup','HybridBackend'});
ValidationSummary=table(max(rel_gap),max(hybrid_gap),min(min_ratio),max(pd_gap), ...
    median(fpi_iter),median(speedup),100*mean(hybrid_backend=="duality"), ...
    'VariableNames',{'MaxDualityCVXRelativeGap','MaxHybridCVXRelativeGap', ...
    'WorstMinSINRRatio','MaxPrimalDualGap','MedianFPIIterations', ...
    'MedianMicrobenchmarkSpeedup','HybridDualityBackendPercent'});

disp(report); disp(ValidationSummary);
matfile=fullfile(results_dir,sprintf('comm_solver_validation_%s.mat',mode));
csvfile=fullfile(results_dir,sprintf('comm_solver_validation_%s.csv',mode));
summary_csv=fullfile(results_dir,sprintf('comm_solver_validation_summary_%s.csv',mode));
save(matfile,'report','ValidationSummary','cases','sigma');
writetable(report,csvfile); writetable(ValidationSummary,summary_csv);
fprintf(['Duality/hybrid validation PASS: max FPI-CVX gap %.3g, max hybrid-CVX gap %.3g, ' ...
    'median micro-benchmark speedup %.1fx.\n'],max(rel_gap),max(hybrid_gap),median(speedup));
fprintf('Validation data saved to %s\n',matfile);
end

function s=compute_sinr(H,W,sigma)
K=size(H,2); s=zeros(K,1);
for k=1:K
    y=H(:,k)'*W;
    s(k)=abs(y(k))^2/(sum(abs(y).^2)-abs(y(k))^2+sigma^2);
end
end
