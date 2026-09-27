function [phi,cost,info] = unit_modulus_rcg(costgrad,phi0,options)
%UNIT_MODULUS_RCG Riemannian conjugate-gradient on the complex circle.
%
% The manifold is M={phi: |phi_n|=1}. costgrad(phi) returns a real cost and
% a complex Euclidean gradient g satisfying df ~= real(g'*dphi).
%
% Polak-Ribiere+ momentum with tangent-space vector transport and periodic
% restart is used. An Armijo backtracking line search guarantees descent of
% the surrogate objective. Compared with plain Riemannian gradient descent,
% RCG typically requires fewer expensive surrogate iterations while keeping
% the same unit-modulus feasibility.

if nargin<3; options=struct(); end
maxiter=get_opt(options,'maxiter',80);
tolgrad=get_opt(options,'tolgradnorm',1e-6);
minstep=get_opt(options,'minstepsize',1e-10);
verbosity=get_opt(options,'verbosity',0);
restart_period=max(1,round(get_opt(options,'restartperiod',12)));
armijo=get_opt(options,'armijo',1e-4);
backtrack=get_opt(options,'backtrack',0.5);
stepgrow=get_opt(options,'stepgrow',1.5);
maxstep=get_opt(options,'maxstep',1.0);

phi=normalize_phase(phi0);
[cost,g]=costgrad(phi); cost=real(cost);
rg=project_tangent(phi,g);
d=-rg;
step_hint=1;

info.cost=nan(maxiter+1,1); info.gradnorm=nan(maxiter+1,1);
info.stepsize=nan(maxiter,1); info.beta=nan(maxiter,1);
info.cost(1)=cost;
last_it=0; accepted_steps=0;

for it=1:maxiter
    last_it=it;
    ng=norm(rg); info.gradnorm(it)=ng;
    if ~isfinite(cost) || ~all(isfinite(rg)) || ng<tolgrad; break; end

    % Ensure a genuine descent direction after numerical transport/momentum.
    slope=real(rg'*d);
    if ~isfinite(slope) || slope>=-1e-12*max(ng*norm(d),1)
        d=-rg; slope=-ng^2;
    end

    step=min(max(step_hint,minstep),maxstep); accepted=false;
    while step>=minstep
        cand=normalize_phase(phi+step*d);
        [cnew,gnew]=costgrad(cand); cnew=real(cnew);
        if isfinite(cnew) && cnew<=cost+armijo*step*slope
            accepted=true; break;
        end
        step=step*backtrack;
    end
    if ~accepted; break; end

    rg_new=project_tangent(cand,gnew);
    rg_old_t=project_tangent(cand,rg);
    d_old_t=project_tangent(cand,d);
    y=rg_new-rg_old_t;
    beta_pr=real(rg_new'*y)/max(real(rg'*rg),1e-18);
    beta=max(beta_pr,0); % PR+ safeguard.
    if mod(it,restart_period)==0 || ~isfinite(beta); beta=0; end
    d_new=-rg_new+beta*d_old_t;
    if real(rg_new'*d_new)>=0
        beta=0; d_new=-rg_new;
    end

    phi=cand; cost=cnew; g=gnew; rg=rg_new; d=d_new;
    accepted_steps=accepted_steps+1;
    info.cost(it+1)=cost; info.stepsize(accepted_steps)=step; info.beta(accepted_steps)=beta;
    step_hint=min(maxstep,max(minstep,step*stepgrow));
    if verbosity>0
        fprintf('RCG %3d cost=% .6e |grad|=% .3e step=% .2e beta=% .2e\n', ...
            it,cost,ng,step,beta);
    end
end

info.iterations=last_it;
last_cost=find(isfinite(info.cost),1,'last');
if isempty(last_cost); last_cost=1; end
info.cost=info.cost(1:last_cost);
last_grad=find(isfinite(info.gradnorm),1,'last');
if isempty(last_grad); info.gradnorm=[]; else; info.gradnorm=info.gradnorm(1:last_grad); end
if accepted_steps==0
    info.stepsize=[]; info.beta=[];
else
    info.stepsize=info.stepsize(1:accepted_steps);
    info.beta=info.beta(1:accepted_steps);
end
end

function rg=project_tangent(phi,g)
rg=g-real(conj(phi).*g).*phi;
end

function phi=normalize_phase(phi)
phi=phi(:); bad=abs(phi)<1e-12; phi(bad)=1; phi=phi./abs(phi);
end

function v=get_opt(s,name,default)
if isstruct(s)&&isfield(s,name)&&~isempty(s.(name)); v=s.(name); else; v=default; end
end
