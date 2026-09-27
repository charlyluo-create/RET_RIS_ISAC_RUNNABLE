function D = comm_solver_stats_delta(A,B)
%COMM_SOLVER_STATS_DELTA Difference of two cumulative telemetry snapshots.
if nargin<1||isempty(A); A=zero_like(); end
if nargin<2||isempty(B); B=comm_solver_telemetry('snapshot'); end
D=B;
fn=fieldnames(D);
for ii=1:numel(fn)
    f=fn{ii};
    if isnumeric(B.(f)) && isfield(A,f) && isnumeric(A.(f))
        D.(f)=B.(f)-A.(f);
        if abs(D.(f))<1e-12; D.(f)=0; end
    end
end
end
function C=zero_like()
C=comm_solver_telemetry('snapshot'); fn=fieldnames(C);
for ii=1:numel(fn); if isnumeric(C.(fn{ii})); C.(fn{ii})=0; end; end
end
