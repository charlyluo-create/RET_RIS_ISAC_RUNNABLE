function [Wc,p_comm,feasible,meta] = solve_comm_min_power_duality(Hu,gamma,sigma_c,opts)
%SOLVE_COMM_MIN_POWER_DUALITY Fast globally-optimal QoS beamforming solver.
%
% Solves
%   min_W ||W||_F^2
%   s.t. SINR_k >= gamma_k,  k=1,...,K
% for the MISO downlink by classical uplink-downlink duality / KKT fixed-point
% iteration. The beam directions follow the virtual-uplink MMSE form and the
% downlink powers are recovered from the SINR equality system.
%
% The routine is self-validating. It checks the achieved SINRs and the
% primal-dual power gap. The calling hybrid wrapper falls back to CVX if any
% check fails, so numerical speed does not come at the expense of correctness.

if nargin<4 || isempty(opts); opts=struct(); end
max_iter=get_opt(opts,'max_iter',500);
tol=get_opt(opts,'tol',1e-7);
damping=get_opt(opts,'damping',1.0);
max_gap=get_opt(opts,'max_primal_dual_gap',2e-4);
sinr_tol=get_opt(opts,'sinr_tol',5e-4);

[M,K]=size(Hu);
Wc=zeros(M,K); p_comm=Inf; feasible=false;
meta=struct('status','Not started','iterations',0,'rel_fixed_point',Inf, ...
    'primal_dual_gap',Inf,'min_sinr_ratio',0,'spectral_radius',Inf, ...
    'dual_power',NaN,'scale',NaN);
if K==0 || M==0 || sigma_c<0 || ~all(isfinite(Hu(:)))
    meta.status='Invalid input'; return;
end
if isscalar(gamma)
    gamma_vec=repmat(real(gamma),K,1);
else
    gamma_vec=real(gamma(:));
end
if numel(gamma_vec)~=K || any(~isfinite(gamma_vec)) || any(gamma_vec<=0)
    meta.status='Invalid gamma'; return;
end

% Normalize channel/noise jointly. This leaves every SINR and the optimal
% transmit power unchanged while improving conditioning for weak channels.
scale=max(norm(Hu,'fro'),1e-12);
Hs=Hu/scale; ns=sigma_c/scale; meta.scale=scale;

% mu_k denotes gamma_k times the conventional SINR dual multiplier. A useful
% positive initialization follows the single-user solution and avoids a
% near-zero first iterate in highly attenuated V2X channels.
hnorm2=sum(abs(Hs).^2,1).';
mu=gamma_vec./max(hnorm2,1e-12);
mu=max(mu,1e-12);
I=eye(M);
rel=Inf;

for it=1:max_iter
    A=I + Hs*diag(mu)*Hs';
    A=(A+A')/2;
    if rcond(A)<1e-14 || any(~isfinite(A(:)))
        meta.status='Ill-conditioned dual matrix'; meta.iterations=it; return;
    end
    X=A\Hs;
    q=real(sum(conj(Hs).*X,1)).';
    if any(~isfinite(q)) || any(q<=1e-14)
        meta.status='Invalid dual response'; meta.iterations=it; return;
    end
    mu_map=gamma_vec./((1+gamma_vec).*q);
    if any(~isfinite(mu_map)) || any(mu_map<=0)
        meta.status='Invalid fixed-point update'; meta.iterations=it; return;
    end
    rel=max(abs(mu_map-mu)./max(abs(mu),1e-12));
    mu_next=(1-damping)*mu + damping*mu_map;
    mu=max(real(mu_next),1e-14);
    if rel<tol; break; end
end
meta.iterations=it; meta.rel_fixed_point=rel;
if rel>max(100*tol,1e-5)
    meta.status='Fixed point did not converge'; return;
end

% Recompute the KKT beam directions at the converged dual point.
A=I + Hs*diag(mu)*Hs'; A=(A+A')/2;
V=A\Hs;
for k=1:K
    nv=norm(V(:,k));
    if ~isfinite(nv) || nv<=1e-14
        meta.status='Degenerate beam direction'; return;
    end
    V(:,k)=V(:,k)/nv;
end

% Downlink power loading: enforce all SINR constraints with equality for the
% fixed optimal directions. p = (I-F)^{-1}u.
F=zeros(K,K); u=zeros(K,1);
for k=1:K
    gains=abs(Hs(:,k)'*V).^2;
    gkk=real(gains(k));
    if ~isfinite(gkk) || gkk<=1e-16
        meta.status='Vanishing desired gain'; return;
    end
    for j=1:K
        if j~=k; F(k,j)=gamma_vec(k)*real(gains(j))/gkk; end
    end
    u(k)=gamma_vec(k)*ns^2/gkk;
end
rho=max(abs(eig(F))); meta.spectral_radius=real(rho);
if ~isfinite(rho) || rho>=1-1e-10 || rcond(eye(K)-F)<1e-13
    meta.status='SINR power-loading system infeasible/ill-conditioned'; return;
end
p=(eye(K)-F)\u;
if any(~isfinite(p)) || any(p<-1e-10)
    meta.status='Invalid downlink powers'; return;
end
p=max(real(p),0);
Wc=V*diag(sqrt(p));

% Match the usual SOCP phase convention h_k^H w_k real and nonnegative.
for k=1:K
    z=Hs(:,k)'*Wc(:,k);
    if abs(z)>0; Wc(:,k)=Wc(:,k)*exp(-1i*angle(z)); end
end
p_comm=real(norm(Wc,'fro')^2);

% Strict original-channel SINR validation.
sinr=zeros(K,1);
for k=1:K
    y=Hu(:,k)'*Wc;
    den=sum(abs(y).^2)-abs(y(k))^2+sigma_c^2;
    sinr(k)=abs(y(k))^2/max(real(den),realmin);
end
ratio=sinr./gamma_vec; meta.min_sinr_ratio=min(ratio);
if any(~isfinite(sinr)) || any(ratio<1-sinr_tol)
    meta.status='SINR validation failed'; Wc=zeros(M,K); p_comm=Inf; return;
end

% Strong-duality certificate for the normalized problem. With mu_k defined
% above, the dual objective is n_s^2 sum_k mu_k.
dual_power=ns^2*sum(mu); meta.dual_power=real(dual_power);
gap=abs(p_comm-dual_power)/max(p_comm,1e-12); meta.primal_dual_gap=real(gap);
if ~isfinite(gap) || gap>max_gap
    meta.status='Primal-dual gap validation failed'; Wc=zeros(M,K); p_comm=Inf; return;
end

feasible=true;
meta.status='Solved by uplink-downlink duality fixed point';
end

function v=get_opt(s,name,default)
if isstruct(s)&&isfield(s,name)&&~isempty(s.(name)); v=s.(name); else; v=default; end
end
