function out = evaluate_external_beamformers(cfg,scen,theta_deg,phi,use_ris,Wc,Wr,gamma)
%EVALUATE_EXTERNAL_BEAMFORMERS Evaluate estimated-CSI beams on true channels.
% Hardware variables and active beamformers are held fixed. The function
% reports actual user SINRs, QoS satisfaction, power, and weighted radar SNR.

out=struct('valid',false,'qos_satisfied',false,'sinr',[],'min_sinr_dB',NaN, ...
    'p_comm',NaN,'p_radar',NaN,'p_total',NaN,'obj_snr',NaN,'theta_deg',NaN);
if nargin<9 || isempty(gamma); gamma=0; end
if isempty(Wc) || size(Wc,1)~=scen.M || size(Wc,2)~=scen.K
    return;
end
if isempty(Wr); Wr=zeros(scen.M,1); end
if size(Wr,1)~=scen.M; return; end

[~,ti]=min(abs(scen.theta_range-theta_deg));
[Hd_u,Hd_t,G]=apply_ret_channels(ti,scen.Hd_user,scen.Hd_tar,scen.G, ...
    scen.gain_u,scen.gain_t,scen.gain_ris);
if use_ris
    phi=normalize_phase(phi,scen.N);
    Hu=build_effective_user_channel(phi,true,Hd_u,scen.Hr_user,G);
else
    phi=[];
    Hu=Hd_u;
end
C=build_radar_matrix(use_ris,phi,Hd_t,scen.Hr_tar,G, ...
    cfg.sigma_r2,cfg.target_weights);

sinr=zeros(scen.K,1);
for k=1:scen.K
    y=Hu(:,k)'*Wc;
    desired=abs(y(k))^2;
    sensing_interference=sum(abs(Hu(:,k)'*Wr).^2);
    interf=sum(abs(y).^2)-desired+sensing_interference;
    sinr(k)=desired/max(interf+cfg.sigma_c2,realmin);
end
W=[Wc,Wr];
out.valid=all(isfinite(sinr)) && all(isfinite(W(:)));
out.sinr=sinr;
out.min_sinr_dB=10*log10(max(min(sinr),realmin));
qos_tol=5e-4;
if isfield(cfg,'csi_qos_tolerance_exp9') && isfinite(cfg.csi_qos_tolerance_exp9)
    qos_tol=cfg.csi_qos_tolerance_exp9;
end
out.qos_satisfied=out.valid && all(sinr>=gamma*(1-qos_tol));
out.p_comm=real(norm(Wc,'fro')^2);
out.p_radar=real(norm(Wr,'fro')^2);
out.p_total=out.p_comm+out.p_radar;
out.obj_snr=real(trace(W'*C*W));
out.theta_deg=scen.theta_range(ti);
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
