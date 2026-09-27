function [mu,lo,hi,n] = masked_mean_ci_linear(X,mask)
%MASKED_MEAN_CI_LINEAR Mean and two-sided Student-t 95% CI under a mask.
% X may be ITER-by-P or ITER-by-P-by-S. MASK may have the same size as X
% or omit the final scheme dimension (ITER-by-P), in which case it is shared
% by every scheme. Non-finite values are excluded after applying MASK.

if nargin<2 || isempty(mask)
    mask=true(size(X));
end
if ndims(X)==2
    X=reshape(X,size(X,1),size(X,2),1);
    was2d=true;
else
    was2d=false;
end
[I,P,S]=size(X);
if isequal(size(mask),[I P])
    mask=repmat(reshape(mask,[I P 1]),[1 1 S]);
elseif ~isequal(size(mask),size(X))
    error('MASK must be ITER-by-P or the same size as X.');
end
mu=nan(P,S); lo=nan(P,S); hi=nan(P,S); n=zeros(P,S);
for p=1:P
    for s=1:S
        v=X(:,p,s);
        ok=mask(:,p,s) & isfinite(v);
        v=v(ok); n(p,s)=numel(v);
        if isempty(v); continue; end
        mu(p,s)=mean(v);
        if numel(v)>1
            half=student_t_critical(0.95,numel(v)-1)*std(v,0)/sqrt(numel(v));
        else
            half=0;
        end
        lo(p,s)=mu(p,s)-half;
        hi(p,s)=mu(p,s)+half;
    end
end
if was2d
    mu=mu(:,1); lo=lo(:,1); hi=hi(:,1); n=n(:,1);
end
end
