function [qhat,lo,hi,n] = bootstrap_quantile_grid(X,q,B,seed)
%BOOTSTRAP_QUANTILE_GRID Bootstrap percentile CI for a grid over dim 1.
% X is ITER-by-P or ITER-by-P-by-Q. Each grid cell is resampled independently.
if nargin<3 || isempty(B); B=1000; end
if nargin<4 || isempty(seed); seed=1; end
if ndims(X)==2
    X=reshape(X,size(X,1),size(X,2),1); was2d=true;
else
    was2d=false;
end
[~,P,Q]=size(X);
qhat=nan(P,Q); lo=nan(P,Q); hi=nan(P,Q); n=zeros(P,Q);
old=rng; cleanup=onCleanup(@()rng(old)); %#ok<NASGU>
rng(seed,'twister');
for p=1:P
    for j=1:Q
        v=X(:,p,j); v=v(isfinite(v)); n(p,j)=numel(v);
        if isempty(v); continue; end
        qhat(p,j)=empirical_quantile(v,q);
        if numel(v)<2 || B<2
            lo(p,j)=qhat(p,j); hi(p,j)=qhat(p,j); continue;
        end
        boot=nan(B,1); nv=numel(v);
        for b=1:B
            ids=randi(nv,nv,1);
            boot(b)=empirical_quantile(v(ids),q);
        end
        lo(p,j)=empirical_quantile(boot,0.025);
        hi(p,j)=empirical_quantile(boot,0.975);
    end
end
if was2d
    qhat=qhat(:,1); lo=lo(:,1); hi=hi(:,1); n=n(:,1);
end
end
