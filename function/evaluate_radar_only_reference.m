function obj = evaluate_radar_only_reference(cfg,scen,theta_deg,phi,P_total,use_ris)
%EVALUATE_RADAR_ONLY_REFERENCE Same-configuration radar-only reference.
% This is conditioned on the supplied RET/RIS configuration and is not a
% global upper bound over passive configurations.
if nargin<6 || isempty(use_ris); use_ris=true; end
[~,ti]=min(abs(scen.theta_range-theta_deg));
[~,Hd_t,G]=apply_ret_channels(ti,scen.Hd_user,scen.Hd_tar,scen.G, ...
    scen.gain_u,scen.gain_t,scen.gain_ris);
if use_ris; phi=normalize_phase(phi,scen.N); else; phi=[]; end
C=build_radar_matrix(use_ris,phi,Hd_t,scen.Hr_tar,G, ...
    cfg.sigma_r2,cfg.target_weights);
[V,D]=eig((C+C')/2); [~,idx]=max(real(diag(D)));
v=V(:,idx); v=v/max(norm(v),eps);
w=sqrt(P_total)*v;
obj=real(w'*C*w);
end

function [Hd_u,Hd_t,Gt]=apply_ret_channels(ti,Hd_u0,Hd_t0,G0,gain_u,gain_t,gain_ris)
Hd_u=Hd_u0; Hd_t=Hd_t0;
for k=1:size(Hd_u0,2); Hd_u(:,k)=sqrt(gain_u(ti,k))*Hd_u0(:,k); end
for t=1:size(Hd_t0,2); Hd_t(:,t)=sqrt(gain_t(ti,t))*Hd_t0(:,t); end
Gt=sqrt(gain_ris(ti))*G0;
end
function phi=normalize_phase(phi,N)
if isempty(phi)||numel(phi)~=N; phi=ones(N,1); else; phi=phi(:); end
bad=abs(phi)<1e-12; phi(bad)=1; phi=phi./abs(phi);
end
