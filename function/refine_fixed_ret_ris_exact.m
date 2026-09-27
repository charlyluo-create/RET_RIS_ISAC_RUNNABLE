function [phi_best,res_best,meta] = refine_fixed_ret_ris_exact( ...
    cfg,scen,theta_deg,phi_seed,res_seed,gamma,P_total)
%REFINE_FIXED_RET_RIS_EXACT Strengthen fixed-RET RIS sensing baselines fairly.
%
% This helper applies the same surrogate-screened / exact-acceptance RIS
% neighborhood idea used by the enhanced joint method, but it never changes
% RET.  It is enabled only for full sensing-oriented baseline optimization;
% the mobility experiment's deliberately restricted fresh-RIS update does not
% call this refinement.

phi_best=normalize_phase(phi_seed,scen.N);
res_best=res_seed;
meta=struct('screened',0,'exact_evals',0,'accepts',0);
if ~isstruct(res_best) || ~isfield(res_best,'feasible') || ~res_best.feasible
    return;
end

[~,ti]=min(abs(scen.theta_range-theta_deg));
[Hd_u,Hd_t,G]=apply_ret_channels(ti,scen.Hd_user,scen.Hd_tar,scen.G, ...
    scen.gain_u,scen.gain_t,scen.gain_ris);
Bsch=get_cfg(cfg,'proposed_polish_blocks_schedule',8);
B=max(1,min(scen.N,round(Bsch(1))));
stepSch=get_cfg(cfg,'proposed_polish_phase_step_schedule',{[22.5 11.25]});
if iscell(stepSch); steps=abs(stepSch{1}); else; steps=abs(stepSch); end
steps=steps(isfinite(steps)&steps>0); if isempty(steps); steps=11.25; end
qsch=get_cfg(cfg,'proposed_polish_exact_topk_max_schedule',2);
topk=max(1,round(qsch(1)));
sweeps=1;
block_id=ceil((1:scen.N)*B/scen.N);

for sw=1:sweeps
    cand_phi=cell(0,1); cand_score=[];
    for bb=1:B
        ids=find(block_id==bb);
        for ss=1:numel(steps)
            for sg=[-1 1]
                q=phi_best;
                q(ids)=q(ids)*exp(1i*sg*steps(ss)*pi/180);
                q=normalize_phase(q,scen.N);
                cand_phi{end+1}=q; %#ok<AGROW>
                cand_score(end+1)=proxy_score(q,Hd_u,Hd_t,G,res_best,cfg,scen,gamma,P_total); %#ok<AGROW>
                meta.screened=meta.screened+1;
            end
        end
    end
    [~,ord]=sort(cand_score,'descend');
    neval=min(topk,sum(isfinite(cand_score))); used=0;
    local_best=res_best; local_phi=phi_best;
    for jj=1:numel(ord)
        if used>=neval; break; end
        qid=ord(jj); if ~isfinite(cand_score(qid)); continue; end
        rr=evaluate_fixed_ret_ris_config(cfg,scen,theta_deg,cand_phi{qid},gamma,P_total,true);
        meta.exact_evals=meta.exact_evals+1; used=used+1;
        if exact_better(rr,local_best,cfg)
            local_best=rr; local_phi=cand_phi{qid};
        end
    end
    if exact_better(local_best,res_best,cfg)
        res_best=local_best; phi_best=local_phi; meta.accepts=meta.accepts+1;
    else
        break;
    end
end
res_best.theta_deg=theta_deg; res_best.phi=phi_best;
end

function score=proxy_score(phi,Hd_u,Hd_t,G,current,cfg,scen,gamma,P_total)
% Same power--channel coupled screening proxy as the enhanced joint method.
Hu=build_effective_user_channel(phi,true,Hd_u,scen.Hr_user,G);
C=build_radar_matrix(true,phi,Hd_t,scen.Hr_tar,G,cfg.sigma_r2,cfg.target_weights);
Z=null(Hu');
if isempty(Z); Geff=0; else; B=Z'*C*Z; B=(B+B')/2; Geff=max(max(real(eig(B))),0); end
a_req=uniform_comm_power_scale(Hu,current.Wc,gamma,sqrt(cfg.sigma_c2));
if ~isfinite(a_req); score=-Inf; return; end
Pbase=real(norm(current.Wc,'fro')^2);
Pcomm_hat=a_req*Pbase;
if ~isfinite(Pcomm_hat) || Pcomm_hat>P_total*(1+1e-9); score=-Inf; return; end
Jcomm_fixed=real(trace(current.Wc'*C*current.Wc));
Jcomm_hat=max(a_req,0)*max(Jcomm_fixed,0);
Pr_hat=max(P_total-Pcomm_hat,0);
score=max(real(Jcomm_hat+Pr_hat*Geff),0);
if ~isfinite(score); score=-Inf; end
end

function a_req=uniform_comm_power_scale(Hu,Wc,gamma,sigma_c)
K=size(Hu,2); noise=sigma_c^2; a_req=0;
if isscalar(gamma); gamma_vec=repmat(gamma,K,1); else; gamma_vec=gamma(:); end
if numel(gamma_vec)~=K; a_req=Inf; return; end
for k=1:K
    hk=Hu(:,k); sig=abs(hk'*Wc(:,k))^2; interf=0;
    for j=1:size(Wc,2); if j~=k; interf=interf+abs(hk'*Wc(:,j))^2; end; end
    margin=sig-gamma_vec(k)*interf;
    if ~isfinite(margin) || margin<=1e-18; a_req=Inf; return; end
    ak=gamma_vec(k)*noise/margin;
    if ~isfinite(ak) || ak<0; a_req=Inf; return; end
    a_req=max(a_req,ak);
end
end

function tf=exact_better(a,b,cfg)
tf=false;
if ~isstruct(a)||~isfield(a,'feasible')||~a.feasible||~isfinite(a.obj_snr); return; end
if ~isstruct(b)||~isfield(b,'feasible')||~b.feasible; tf=true; return; end
ratio=max(1+get_cfg(cfg,'proposed_refine_accept_tol',1e-10), ...
    10^(get_cfg(cfg,'proposed_polish_min_gain_dB',0)/10));
tf=a.obj_snr>b.obj_snr*ratio;
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

function v=get_cfg(cfg,name,default)
if isfield(cfg,name)&&~isempty(cfg.(name)); v=cfg.(name); else; v=default; end
end
