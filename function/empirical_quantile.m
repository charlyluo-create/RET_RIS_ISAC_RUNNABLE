function qv = empirical_quantile(x,p)
%EMPIRICAL_QUANTILE Toolbox-free linearly interpolated sample quantile.
% X may be any real numeric array; non-finite entries are removed. P must
% contain probabilities in [0,1]. The interpolation index is 1+(n-1)P.
if nargin<2; error('Both x and p are required.'); end
if ~isnumeric(x) || ~isnumeric(p); error('x and p must be numeric.'); end
p_shape=size(p); p=p(:);
if any(~isfinite(p)|p<0|p>1); error('p must lie in [0,1].'); end
x=x(:); x=x(isfinite(x)); x=sort(x);
qv=nan(size(p)); n=numel(x);
if n==0; qv=reshape(qv,p_shape); return; end
if n==1; qv(:)=x; qv=reshape(qv,p_shape); return; end
for ii=1:numel(p)
    h=1+(n-1)*p(ii); lo=floor(h); hi=ceil(h);
    if lo==hi; qv(ii)=x(lo);
    else; qv(ii)=x(lo)+(h-lo)*(x(hi)-x(lo)); end
end
qv=reshape(qv,p_shape);
end
