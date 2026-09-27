function v = comm_solver_stats_summary(S)
%COMM_SOLVER_STATS_SUMMARY Convert telemetry struct to a compact numeric row.
% Columns: calls, duality successes, CVX attempts, hybrid fallbacks,
% fallback rate, mean FPI iterations, mean primal-dual gap, solver wall time.
if nargin<1||isempty(S)
    S=comm_solver_telemetry('snapshot'); fn=fieldnames(S);
    for ii=1:numel(fn); if isnumeric(S.(fn{ii})); S.(fn{ii})=0; end; end
end
calls=max(S.total_calls,0);
fbRate=S.hybrid_fallbacks/max(S.requested_hybrid,1);
if S.fpi_iterations_count>0; meanIt=S.fpi_iterations_sum/S.fpi_iterations_count; else; meanIt=NaN; end
if S.primal_dual_gap_count>0; meanGap=S.primal_dual_gap_sum/S.primal_dual_gap_count; else; meanGap=NaN; end
v=[calls,S.duality_successes,S.cvx_attempts,S.hybrid_fallbacks, ...
   fbRate,meanIt,meanGap,S.total_wall_s];
end
