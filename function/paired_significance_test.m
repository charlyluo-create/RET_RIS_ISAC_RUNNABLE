function stats = paired_significance_test(A,B,mode,alpha)
%PAIRED_SIGNIFICANCE_TEST Toolbox-free paired effect and significance tests.
% MODE='db_ratio' uses 10log10(A/B); MODE='difference' uses A-B.
% Reports a paired t-test and a Wilcoxon signed-rank normal approximation.

if nargin<3 || isempty(mode); mode='difference'; end
if nargin<4 || isempty(alpha); alpha=0.05; end
A=A(:); B=B(:);
ok=isfinite(A)&isfinite(B);
switch lower(mode)
    case 'db_ratio'
        ok=ok&A>0&B>0;
        d=10*log10(A(ok)./B(ok));
    case 'difference'
        d=A(ok)-B(ok);
    otherwise
        error('Unknown mode: %s',mode);
end
n=numel(d);
stats=struct('n',n,'mean',NaN,'median',NaN,'sd',NaN,'ci_low',NaN, ...
    'ci_high',NaN,'t_stat',NaN,'p_ttest',NaN,'wilcoxon_z',NaN, ...
    'p_wilcoxon',NaN,'alpha',alpha,'significant_t',false, ...
    'significant_wilcoxon',false,'samples',d);
if n==0; return; end
stats.mean=mean(d); stats.median=median(d); stats.sd=std(d,0);
if n==1
    stats.ci_low=d; stats.ci_high=d; return;
end
se=stats.sd/sqrt(n);
tcrit=student_t_critical(1-alpha,n-1);
stats.ci_low=stats.mean-tcrit*se;
stats.ci_high=stats.mean+tcrit*se;
if se>0
    t=stats.mean/se; nu=n-1;
    p=betainc(nu/(nu+t^2),nu/2,0.5);
    stats.t_stat=t; stats.p_ttest=min(max(p,0),1);
else
    stats.t_stat=sign(stats.mean)*Inf;
    stats.p_ttest=double(abs(stats.mean)<eps);
end
stats.significant_t=stats.p_ttest<alpha;

% Wilcoxon signed-rank test with zero differences removed and tie correction.
zv=d(abs(d)>1e-14);
nw=numel(zv);
if nw>=2
    av=abs(zv); r=tied_ranks_local(av);
    Wp=sum(r(zv>0));
    mu=nw*(nw+1)/4;
    [~,~,ic]=unique(av);
    counts=accumarray(ic,1);
    counts=counts(counts>1);
    tie_term=sum(counts.*(counts+1).*(2*counts+1));
    varW=nw*(nw+1)*(2*nw+1)/24-tie_term/48;
    if varW>0
        corr=0.5*sign(Wp-mu);
        z=(Wp-mu-corr)/sqrt(varW);
        pw=erfc(abs(z)/sqrt(2));
        stats.wilcoxon_z=z;
        stats.p_wilcoxon=min(max(pw,0),1);
        stats.significant_wilcoxon=stats.p_wilcoxon<alpha;
    end
end
end

function r=tied_ranks_local(x)
[xs,ord]=sort(x(:)); n=numel(xs); rs=zeros(n,1); i=1;
while i<=n
    j=i;
    while j<n && abs(xs(j+1)-xs(i))<=1e-12*max(1,abs(xs(i))); j=j+1; end
    rs(i:j)=0.5*(i+j); i=j+1;
end
r=zeros(n,1); r(ord)=rs;
end
