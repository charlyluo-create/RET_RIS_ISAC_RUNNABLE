function [mu,lo,hi,n,samples] = paired_db_statistics(A,B,mask)
%PAIRED_DB_STATISTICS Statistics of 10log10(A/B) on paired realizations.
% A and B must have equal size (ITER-by-P or ITER-by-P-by-S). MASK is
% optional and may be shared across the final dimension.
if ~isequal(size(A),size(B)); error('A and B must have equal size.'); end
if nargin<3 || isempty(mask); mask=true(size(A)); end
if ndims(A)==2
    A=reshape(A,size(A,1),size(A,2),1);
    B=reshape(B,size(B,1),size(B,2),1);
    was2d=true;
else
    was2d=false;
end
[I,P,S]=size(A);
if isequal(size(mask),[I P])
    mask=repmat(reshape(mask,[I P 1]),[1 1 S]);
elseif ~isequal(size(mask),size(A))
    error('MASK must be ITER-by-P or the same size as A/B.');
end
mu=nan(P,S); lo=nan(P,S); hi=nan(P,S); n=zeros(P,S);
samples=nan(I,P,S);
for p=1:P
    for s=1:S
        ok=mask(:,p,s) & isfinite(A(:,p,s)) & isfinite(B(:,p,s)) & ...
            A(:,p,s)>0 & B(:,p,s)>0;
        d=10*log10(A(ok,p,s)./B(ok,p,s));
        samples(1:numel(d),p,s)=d;
        n(p,s)=numel(d);
        if isempty(d); continue; end
        mu(p,s)=mean(d);
        if numel(d)>1; half=student_t_critical(0.95,numel(d)-1)*std(d,0)/sqrt(numel(d)); else; half=0; end
        lo(p,s)=mu(p,s)-half; hi(p,s)=mu(p,s)+half;
    end
end
if was2d
    mu=mu(:,1); lo=lo(:,1); hi=hi(:,1); n=n(:,1); samples=samples(:,:,1);
end
end
