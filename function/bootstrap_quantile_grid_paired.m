function [qhat,lo,hi,n,bootq] = bootstrap_quantile_grid_paired(X,q,B,seed)
%BOOTSTRAP_QUANTILE_GRID_PAIRED Paired bootstrap percentile CI over dim 1.
%
% Unlike bootstrap_quantile_grid, one Monte-Carlo index vector is drawn per
% bootstrap replicate and reused for every grid cell. This preserves the
% common-random-number pairing across speeds/update periods and therefore
% supports paired comparisons of neighboring operating points.
%
% X     : ITER-by-P or ITER-by-P-by-Q array.
% q     : requested quantile in [0,1].
% B     : number of bootstrap replicates.
% seed  : deterministic bootstrap seed.
% bootq : B-by-P-by-Q bootstrap quantile surfaces.

if nargin<3 || isempty(B); B=1000; end
if nargin<4 || isempty(seed); seed=1; end
validateattributes(q,{'numeric'},{'scalar','>=',0,'<=',1});
validateattributes(B,{'numeric'},{'scalar','integer','>=',1});

if ndims(X)==2
    X=reshape(X,size(X,1),size(X,2),1); was2d=true;
else
    was2d=false;
end
[I,P,Q]=size(X);
qhat=nan(P,Q); lo=nan(P,Q); hi=nan(P,Q); n=zeros(P,Q);
for p=1:P
    for j=1:Q
        v=X(:,p,j); v=v(isfinite(v)); n(p,j)=numel(v);
        if ~isempty(v); qhat(p,j)=empirical_quantile(v,q); end
    end
end

bootq=nan(B,P,Q);
old=rng; cleanup=onCleanup(@()rng(old)); %#ok<NASGU>
rng(seed,'twister');
for b=1:B
    % The same resampled MC identities are used across the full surface.
    ids=randi(I,I,1);
    Xb=X(ids,:,:);
    for p=1:P
        for j=1:Q
            v=Xb(:,p,j); v=v(isfinite(v));
            if ~isempty(v); bootq(b,p,j)=empirical_quantile(v,q); end
        end
    end
end
for p=1:P
    for j=1:Q
        v=bootq(:,p,j); v=v(isfinite(v));
        if isempty(v); continue; end
        lo(p,j)=empirical_quantile(v,0.025);
        hi(p,j)=empirical_quantile(v,0.975);
    end
end

if was2d
    qhat=qhat(:,1); lo=lo(:,1); hi=hi(:,1); n=n(:,1); bootq=bootq(:,:,1);
end
end
