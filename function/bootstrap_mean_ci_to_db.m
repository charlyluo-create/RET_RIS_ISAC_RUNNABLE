function [mu_db,lo_db,hi_db,n,mu_linear] = bootstrap_mean_ci_to_db( ...
    X,mask,reference,B,seed)
%BOOTSTRAP_MEAN_CI_TO_DB Percentile-bootstrap CI for a nonnegative mean.
% Means are computed in the linear domain and then displayed in dB/dBm.
if nargin<3 || isempty(reference); reference=1; end
if nargin<4 || isempty(B); B=1000; end
if nargin<5 || isempty(seed); seed=1; end
if ndims(X)==2
    X=reshape(X,size(X,1),size(X,2),1); was2d=true;
else
    was2d=false;
end
[I,P,S]=size(X);
if isequal(size(mask),[I P])
    mask=repmat(reshape(mask,[I P 1]),[1 1 S]);
elseif ~isequal(size(mask),size(X))
    error('MASK must be ITER-by-P or the same size as X.');
end
mu_linear=nan(P,S); lo_linear=nan(P,S); hi_linear=nan(P,S); n=zeros(P,S);
state=rng; cleanup=onCleanup(@()rng(state)); %#ok<NASGU>
rng(seed,'twister');
for p=1:P
    for s=1:S
        v=X(:,p,s); ok=mask(:,p,s)&isfinite(v)&v>=0; v=v(ok);
        n(p,s)=numel(v); if isempty(v); continue; end
        mu_linear(p,s)=mean(v);
        if numel(v)==1 || B<=1
            lo_linear(p,s)=mu_linear(p,s); hi_linear(p,s)=mu_linear(p,s);
            continue;
        end
        bm=zeros(B,1);
        nv=numel(v);
        for b=1:B
            idx=randi(nv,nv,1); bm(b)=mean(v(idx));
        end
        lo_linear(p,s)=empirical_quantile(bm,0.025);
        hi_linear(p,s)=empirical_quantile(bm,0.975);
    end
end
mu_db=to_db(mu_linear,reference);
lo_db=to_db(lo_linear,reference);
hi_db=to_db(hi_linear,reference);
if was2d
    mu_db=mu_db(:,1); lo_db=lo_db(:,1); hi_db=hi_db(:,1);
    n=n(:,1); mu_linear=mu_linear(:,1);
end
end

function y=to_db(x,reference)
y=nan(size(x)); valid=isfinite(x)&x>0;
y(valid)=10*log10(x(valid)/reference);
% A zero lower percentile is represented as NaN rather than an artificial
% thousands-of-dB error bar. QoS rate is reported beside this metric.
end
