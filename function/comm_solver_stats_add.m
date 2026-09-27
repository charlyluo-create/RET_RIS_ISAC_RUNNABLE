function C = comm_solver_stats_add(varargin)
%COMM_SOLVER_STATS_ADD Add disjoint telemetry blocks without mutating telemetry.
C=zero_like();
for kk=1:nargin
    A=varargin{kk}; if isempty(A)||~isstruct(A); continue; end
    fn=fieldnames(C);
    for ii=1:numel(fn)
        f=fn{ii};
        if isfield(A,f)&&isnumeric(A.(f)); C.(f)=C.(f)+A.(f); end
    end
end
end
function C=zero_like()
S=comm_solver_telemetry('snapshot'); C=S; fn=fieldnames(C);
for ii=1:numel(fn); if isnumeric(C.(fn{ii})); C.(fn{ii})=0; end; end
end
