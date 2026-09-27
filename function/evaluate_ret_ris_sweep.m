function sweep = evaluate_ret_ris_sweep(cfg,scen,gamma_values,P_total_values,solver_seed,opts)
%EVALUATE_RET_RIS_SWEEP Run the unified five-scheme evaluator over one axis.
%
% This is the single execution path used by the total-power sweep (Exp. 2)
% and QoS sweep (Exp. 3).  GAMMA_VALUES and P_TOTAL_VALUES may be scalars or
% equal-length vectors.  Continuation uses the preceding point only as a
% candidate warm start; every accepted solution is still selected by the
% exact MinPc or weighted radar-SNR objective in EVALUATE_RET_RIS_SCHEMES.
%
% Optional OPTS fields:
%   use_continuation (default true)
%   point_seed_stride (used by seed_policy='sequential')
%   seed_policy ('operating_point' by default, or 'sequential')
%   evaluator_opts (additional options forwarded to the common evaluator)

if nargin<6 || isempty(opts); opts=struct(); end
if nargin<5 || isempty(solver_seed); solver_seed=cfg.seed_base; end
use_continuation=get_opt(opts,'use_continuation',true);
point_seed_stride=get_opt(opts,'point_seed_stride', ...
    get_cfg(cfg,'solver_point_seed_stride',104729));
seed_policy=lower(char(get_opt(opts,'seed_policy','operating_point')));
base_opts=get_opt(opts,'evaluator_opts',struct());

gamma_values=gamma_values(:).';
P_total_values=P_total_values(:).';
L=max(numel(gamma_values),numel(P_total_values));
if numel(gamma_values)==1; gamma_values=repmat(gamma_values,1,L); end
if numel(P_total_values)==1; P_total_values=repmat(P_total_values,1,L); end
if numel(gamma_values)~=L || numel(P_total_values)~=L
    error('gamma_values and P_total_values must be scalar or have equal length.');
end

S=5;
sweep.Obj=nan(L,S); sweep.Pcomm=nan(L,S); sweep.Pradar=nan(L,S);
sweep.Feasible=false(L,S); sweep.Theta=nan(L,S);
sweep.G_eff=nan(L,S); sweep.Jcomm=nan(L,S); sweep.Jradar=nan(L,S);
sweep.Phi=cell(L,S); sweep.History=cell(L,S);
sweep.RuntimeBlock=nan(L,S); sweep.RuntimeCumulative=nan(L,S);
sweep.gamma=gamma_values; sweep.P_total=P_total_values;
sweep.solver_seed=zeros(L,1);

phi_fixed=[]; phi_minpc=[]; phi_radar=[];
theta_minpc=[]; theta_radar=[];
for q=1:L
    eopts=base_opts;
    if use_continuation && q>1
        eopts.phi_seed_fixed_minpc=phi_fixed;
        eopts.phi_seed_minpc=phi_minpc;
        eopts.phi_seed_radar=phi_radar;
        eopts.theta_seed_minpc_deg=theta_minpc;
        eopts.theta_seed_radar_deg=theta_radar;
    end
    switch seed_policy
        case 'operating_point'
            this_seed=ret_ris_operating_point_seed(solver_seed,gamma_values(q),P_total_values(q));
        case 'sequential'
            this_seed=double(solver_seed)+(q-1)*double(point_seed_stride);
        otherwise
            error('Unknown seed_policy: %s',seed_policy);
    end
    sweep.solver_seed(q)=this_seed;
    rng(this_seed,'twister');
    r=evaluate_ret_ris_schemes(cfg,scen,gamma_values(q),P_total_values(q),eopts);
    sweep.Obj(q,:)=r.Obj; sweep.Pcomm(q,:)=r.Pcomm;
    sweep.Pradar(q,:)=r.Pradar; sweep.Feasible(q,:)=r.Feasible;
    sweep.Theta(q,:)=r.Theta; sweep.G_eff(q,:)=r.G_eff;
    sweep.Jcomm(q,:)=r.Jcomm; sweep.Jradar(q,:)=r.Jradar;
    sweep.Phi(q,:)=r.Phi; sweep.History(q,:)=r.History;
    sweep.RuntimeBlock(q,:)=r.RuntimeBlock;
    sweep.RuntimeCumulative(q,:)=r.RuntimeCumulative;

    if use_continuation
        if r.Feasible(2) && ~isempty(r.Phi{2}); phi_fixed=r.Phi{2}; end
        if r.Feasible(4)
            if ~isempty(r.Phi{4}); phi_minpc=r.Phi{4}; end
            theta_minpc=r.Theta(4);
        end
        if r.Feasible(5)
            if ~isempty(r.Phi{5}); phi_radar=r.Phi{5}; end
            theta_radar=r.Theta(5);
        end
    end
end
end

function value=get_opt(s,name,default)
if isstruct(s) && isfield(s,name) && ~isempty(s.(name)); value=s.(name);
else; value=default; end
end

function value=get_cfg(s,name,default)
if isstruct(s) && isfield(s,name) && ~isempty(s.(name)); value=s.(name);
else; value=default; end
end
