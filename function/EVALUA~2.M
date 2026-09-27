function res = evaluate_fixed_split_configuration(cfg,scen,theta_deg,phi, ...
    gamma,P_total,alpha,use_ris)
%EVALUATE_FIXED_SPLIT_CONFIGURATION Fair full-budget power-split ablation.
% The passive RET/RIS state is fixed.  The minimum-power communication beams
% are uniformly scaled to use (1-alpha)P_total and the dedicated sensing beam
% uses alpha P_total in Null(Hu').  Hence no total-power budget is discarded.

if nargin<8 || isempty(use_ris); use_ris=true; end
res=empty_result();
if ~isscalar(alpha) || alpha<0 || alpha>1; error('alpha must lie in [0,1].'); end
base=evaluate_fixed_ret_ris_config(cfg,scen,theta_deg,phi,gamma,P_total,use_ris);
if ~base.feasible; return; end
P_comm_budget=(1-alpha)*P_total;
P_radar=alpha*P_total;
if base.p_comm>P_comm_budget*(1+1e-6); return; end
if base.p_comm>1e-14
    Wc=sqrt(P_comm_budget/base.p_comm)*base.Wc;
else
    Wc=base.Wc;
end
Wr=radar_beam_with_fixed_power(base.Hu,base.C,P_radar);
W=[Wc,Wr];
res.feasible=true; res.p_comm=P_comm_budget; res.p_radar=P_radar;
res.obj_snr=real(trace(W'*base.C*W));
res.Wc=Wc; res.Wr=Wr; res.Hu=base.Hu; res.C=base.C;
res.theta_deg=base.theta_deg;
if use_ris; res.phi=normalize_phase(phi,scen.N); else; res.phi=[]; end
end

function Wr=radar_beam_with_fixed_power(Hu,C,P_radar)
M=size(Hu,1); Wr=zeros(M,1);
if P_radar<=1e-12; return; end
Z=null(Hu'); if isempty(Z); return; end
B=(Z'*C*Z); B=(B+B')/2;
[V,D]=eig(B); [~,idx]=max(real(diag(D)));
v=Z*V(:,idx); if norm(v)<1e-12; v=Z(:,1); end
Wr=sqrt(P_radar)*v/max(norm(v),eps);
end

function phi=normalize_phase(phi,N)
if isempty(phi)||numel(phi)~=N; phi=ones(N,1); else; phi=phi(:); end
bad=abs(phi)<1e-12; phi(bad)=1; phi=phi./abs(phi);
end

function res=empty_result()
res=struct('feasible',false,'p_comm',NaN,'p_radar',NaN,'obj_snr',NaN, ...
    'Wc',[],'Wr',[],'Hu',[],'C',[],'theta_deg',NaN,'phi',[]);
end
