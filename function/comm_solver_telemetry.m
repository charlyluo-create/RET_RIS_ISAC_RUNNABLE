function out = comm_solver_telemetry(action,payload)
%COMM_SOLVER_TELEMETRY Lightweight per-worker instrumentation for QoS solves.
%
% Commands:
%   comm_solver_telemetry('reset')
%   S = comm_solver_telemetry('snapshot')
%   comm_solver_telemetry('record',recordStruct)
%
% The accumulator is intentionally worker-local (persistent state). This is
% useful inside PARFOR because each Monte-Carlo realization can reset, measure
% and return its own solver statistics without synchronization overhead.

persistent S
if isempty(S); S=zero_stats(); end
if nargin<1 || isempty(action); action='snapshot'; end
if nargin<2; payload=struct(); end
switch lower(strtrim(char(action)))
    case 'reset'
        S=zero_stats(); out=S;
    case {'snapshot','get'}
        out=S;
    case 'record'
        r=payload;
        S.total_calls=S.total_calls+1;
        S.total_wall_s=S.total_wall_s+getnum(r,'total_time_s',0);
        S.final_feasible_calls=S.final_feasible_calls+double(getbool(r,'final_feasible',false));
        mode=lower(strtrim(char(getfielddef(r,'requested_mode',''))));
        if strcmp(mode,'hybrid'); S.requested_hybrid=S.requested_hybrid+1; end
        if any(strcmp(mode,{'duality','fast'})); S.requested_duality=S.requested_duality+1; end
        if any(strcmp(mode,{'cvx','socp'})); S.requested_cvx=S.requested_cvx+1; end
        da=getbool(r,'duality_attempted',false); ds=getbool(r,'duality_success',false);
        ca=getbool(r,'cvx_attempted',false); cs=getbool(r,'cvx_success',false);
        fb=getbool(r,'hybrid_fallback',false);
        S.duality_attempts=S.duality_attempts+double(da);
        S.duality_successes=S.duality_successes+double(ds);
        S.duality_failures=S.duality_failures+double(da && ~ds);
        S.cvx_attempts=S.cvx_attempts+double(ca);
        S.cvx_successes=S.cvx_successes+double(cs);
        S.cvx_failures=S.cvx_failures+double(ca && ~cs);
        S.hybrid_fallbacks=S.hybrid_fallbacks+double(fb);
        S.duality_wall_s=S.duality_wall_s+getnum(r,'duality_time_s',0);
        S.cvx_wall_s=S.cvx_wall_s+getnum(r,'cvx_time_s',0);
        it=getnum(r,'fpi_iterations',NaN);
        if isfinite(it)
            S.fpi_iterations_sum=S.fpi_iterations_sum+it;
            S.fpi_iterations_count=S.fpi_iterations_count+1;
        end
        gap=getnum(r,'primal_dual_gap',NaN);
        if isfinite(gap)
            S.primal_dual_gap_sum=S.primal_dual_gap_sum+gap;
            S.primal_dual_gap_count=S.primal_dual_gap_count+1;
        end
        sr=getnum(r,'min_sinr_ratio',NaN);
        if isfinite(sr)
            S.min_sinr_ratio_sum=S.min_sinr_ratio_sum+sr;
            S.min_sinr_ratio_count=S.min_sinr_ratio_count+1;
        end
        out=S;
    otherwise
        error('Unknown comm_solver_telemetry action: %s',action);
end
end

function S=zero_stats()
S=struct( ...
    'total_calls',0,'total_wall_s',0,'final_feasible_calls',0, ...
    'requested_hybrid',0,'requested_duality',0,'requested_cvx',0, ...
    'duality_attempts',0,'duality_successes',0,'duality_failures',0, ...
    'cvx_attempts',0,'cvx_successes',0,'cvx_failures',0, ...
    'hybrid_fallbacks',0,'duality_wall_s',0,'cvx_wall_s',0, ...
    'fpi_iterations_sum',0,'fpi_iterations_count',0, ...
    'primal_dual_gap_sum',0,'primal_dual_gap_count',0, ...
    'min_sinr_ratio_sum',0,'min_sinr_ratio_count',0);
end
function v=getnum(s,name,default)
v=getfielddef(s,name,default); if ~isnumeric(v)||~isscalar(v); v=default; end
end
function v=getbool(s,name,default)
v=getfielddef(s,name,default); if isempty(v); v=default; end; v=logical(v(1));
end
function v=getfielddef(s,name,default)
if isstruct(s)&&isfield(s,name)&&~isempty(s.(name)); v=s.(name); else; v=default; end
end
