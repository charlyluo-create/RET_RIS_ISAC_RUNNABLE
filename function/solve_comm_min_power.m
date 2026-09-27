function [Wc,p_comm,feasible,status_out,meta_out] = solve_comm_min_power(Hu,gamma,sigma_c,opts)
%SOLVE_COMM_MIN_POWER Hybrid fast/exact minimum-power downlink beamforming.
%
% R2 implementation note (MATLAB R2025b compatible):
%   The CVX model is deliberately isolated in solve_comm_min_power_cvx_backend.m.
%   This avoids a MATLAB static-workspace conflict between CVX's dynamic
%   variable creation and nested telemetry callbacks. The numerical problem
%   solved by the CVX fallback is unchanged.
%
% Default mode is 'hybrid':
%   1) solve the QoS problem by uplink-downlink-duality fixed-point iteration;
%   2) validate the original SINRs and primal-dual power gap;
%   3) fall back to the CVX/SOCP backend only if validation fails.
%
% Set RET_RIS_COMM_SOLVER to 'hybrid' (default), 'duality', or 'cvx'.
% An optional fourth struct can override this through opts.mode.
%
% The optional fifth output META_OUT exposes backend/fallback/timing data.
% Every completed call is recorded by comm_solver_telemetry for runtime and
% complexity auditing only; telemetry never changes the optimizer result.

if nargin<4 || isempty(opts); opts=struct(); end
mode=getenv('RET_RIS_COMM_SOLVER');
if isfield(opts,'mode') && ~isempty(opts.mode); mode=char(opts.mode); end
if isempty(mode); mode='hybrid'; end
mode=lower(strtrim(mode));

[M,K]=size(Hu);
Wc=zeros(M,K); p_comm=Inf; feasible=false; status_out='Not started';
meta_out=struct('requested_mode',mode,'backend_used','none', ...
    'duality_attempted',false,'duality_success',false, ...
    'cvx_attempted',false,'cvx_success',false,'hybrid_fallback',false, ...
    'fpi_iterations',NaN,'primal_dual_gap',NaN,'min_sinr_ratio',NaN, ...
    'duality_time_s',0,'cvx_time_s',0,'total_time_s',NaN, ...
    'final_feasible',false,'status','Not started');
t_call=tic;

% ---------- Input validation ----------
if isempty(Hu) || sigma_c<0 || ~all(isfinite(Hu(:)))
    status_out='Invalid input';
    meta_out=finalize_comm_meta(meta_out,status_out,feasible,t_call);
    return;
end
if isscalar(gamma)
    gamma_vec=repmat(gamma,K,1);
else
    gamma_vec=gamma(:);
end
if numel(gamma_vec)~=K || any(gamma_vec<=0) || any(~isfinite(gamma_vec))
    status_out='Invalid gamma';
    meta_out=finalize_comm_meta(meta_out,status_out,feasible,t_call);
    return;
end

% ---------- Fast exact-structure backend: uplink-downlink duality ----------
duality_failure='not attempted';
if any(strcmp(mode,{'hybrid','duality','fast'}))
    meta_out.duality_attempted=true;
    td=tic;
    try
        [Wd,pd,fd,meta_d]=solve_comm_min_power_duality(Hu,gamma_vec,sigma_c,opts);
    catch ME
        Wd=zeros(M,K); pd=Inf; fd=false;
        meta_d=struct('iterations',NaN,'primal_dual_gap',NaN, ...
            'min_sinr_ratio',NaN,'status',['exception: ' ME.message]);
    end
    meta_out.duality_time_s=toc(td);
    if isfield(meta_d,'iterations');       meta_out.fpi_iterations=meta_d.iterations; end
    if isfield(meta_d,'primal_dual_gap'); meta_out.primal_dual_gap=meta_d.primal_dual_gap; end
    if isfield(meta_d,'min_sinr_ratio');  meta_out.min_sinr_ratio=meta_d.min_sinr_ratio; end

    if fd
        Wc=Wd; p_comm=pd; feasible=true;
        meta_out.duality_success=true;
        meta_out.backend_used='duality';
        status_out=sprintf('Duality-FPI solved (it=%g, gap=%.2e)', ...
            meta_out.fpi_iterations,meta_out.primal_dual_gap);
        meta_out=finalize_comm_meta(meta_out,status_out,feasible,t_call);
        return;
    end

    if isfield(meta_d,'status') && ~isempty(meta_d.status)
        duality_failure=char(meta_d.status);
    else
        duality_failure='unspecified duality failure';
    end

    if strcmp(mode,'duality') || strcmp(mode,'fast')
        status_out=['Duality-FPI failed: ' duality_failure];
        meta_out=finalize_comm_meta(meta_out,status_out,feasible,t_call);
        return;
    end
end

% ---------- Solver-mode validation ----------
if ~any(strcmp(mode,{'hybrid','cvx','socp'}))
    status_out=['Unknown solver mode: ' mode];
    meta_out=finalize_comm_meta(meta_out,status_out,feasible,t_call);
    return;
end

% ---------- CVX/SOCP fallback / reproducibility backend ----------
meta_out.cvx_attempted=true;
tc=tic;
[Wc_cvx,p_cvx,fea_cvx,status_cvx,meta_cvx] = ...
    solve_comm_min_power_cvx_backend(Hu,gamma_vec,sigma_c);
meta_out.cvx_time_s=toc(tc);

if isfield(meta_cvx,'min_sinr_ratio')
    meta_out.min_sinr_ratio=meta_cvx.min_sinr_ratio;
end
meta_out.cvx_success=logical(fea_cvx);

if fea_cvx
    Wc=Wc_cvx; p_comm=p_cvx; feasible=true;
    meta_out.backend_used='cvx';
    if strcmp(mode,'hybrid')
        meta_out.hybrid_fallback=true;
        status_out=sprintf('%s (hybrid fallback after: %s)',status_cvx,duality_failure);
    else
        status_out=status_cvx;
    end
else
    if strcmp(mode,'hybrid')
        meta_out.hybrid_fallback=true;
        status_out=sprintf('%s (hybrid fallback after: %s)',status_cvx,duality_failure);
    else
        status_out=status_cvx;
    end
end

meta_out=finalize_comm_meta(meta_out,status_out,feasible,t_call);
end

function meta_out=finalize_comm_meta(meta_out,status_out,feasible,t_call)
% Separate local subfunction (not nested) keeps the parent workspace dynamic,
% which is essential for CVX compatibility on recent MATLAB releases.
meta_out.total_time_s=toc(t_call);
meta_out.final_feasible=logical(feasible);
meta_out.status=status_out;
try
    comm_solver_telemetry('record',meta_out);
catch
    % Telemetry is diagnostic only and must never break optimization.
end
end
