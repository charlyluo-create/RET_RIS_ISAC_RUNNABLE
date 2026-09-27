function S = paired_adjacent_trend_stats(X,qhat,bootq,xvals)
%PAIRED_ADJACENT_TREND_STATS Paired adjacent-period trend diagnostics.
%
% X is ITER-by-V-by-T and must preserve MC identity across T. qhat is V-by-T
% (typically the empirical P90), and bootq is B-by-V-by-T from a paired
% bootstrap. The routine never enforces monotonicity; it only reports whether
% an observed local decrease is statistically resolved by paired 95% CIs.

[I,V,T]=size(X); %#ok<ASGLU>
if numel(xvals)~=T; error('xvals length must equal size(X,3).'); end
S=struct();
S.x_left=xvals(1:end-1); S.x_right=xvals(2:end);
S.MeanDelta=nan(V,T-1); S.MeanDeltaLo=nan(V,T-1); S.MeanDeltaHi=nan(V,T-1);
S.MedianDelta=nan(V,T-1); S.NondecreasePct=nan(V,T-1); S.PairCount=zeros(V,T-1);
S.QuantileDelta=nan(V,T-1); S.QuantileDeltaLo=nan(V,T-1); S.QuantileDeltaHi=nan(V,T-1);
for v=1:V
    for t=1:T-1
        a=X(:,v,t); b=X(:,v,t+1); mask=isfinite(a)&isfinite(b);
        d=b(mask)-a(mask); S.PairCount(v,t)=numel(d);
        if ~isempty(d)
            S.MeanDelta(v,t)=mean(d); S.MedianDelta(v,t)=median(d);
            S.NondecreasePct(v,t)=100*mean(d>=0);
            if numel(d)>1
                hw=student_t_critical(0.95,numel(d)-1)*std(d,0)/sqrt(numel(d));
            else
                hw=0;
            end
            S.MeanDeltaLo(v,t)=S.MeanDelta(v,t)-hw;
            S.MeanDeltaHi(v,t)=S.MeanDelta(v,t)+hw;
        end
        if all(size(qhat)>=[v t+1]) && isfinite(qhat(v,t)) && isfinite(qhat(v,t+1))
            S.QuantileDelta(v,t)=qhat(v,t+1)-qhat(v,t);
        end
        if ~isempty(bootq)
            bd=bootq(:,v,t+1)-bootq(:,v,t); bd=bd(isfinite(bd));
            if ~isempty(bd)
                S.QuantileDeltaLo(v,t)=empirical_quantile(bd,0.025);
                S.QuantileDeltaHi(v,t)=empirical_quantile(bd,0.975);
            end
        end
    end
end
S.SignificantMeanDecrease = S.MeanDeltaHi < 0;
S.SignificantQuantileDecrease = S.QuantileDeltaHi < 0;
S.NumSignificantMeanDecreases = nnz(S.SignificantMeanDecrease);
S.NumSignificantQuantileDecreases = nnz(S.SignificantQuantileDecrease);
S.NumPointMeanDecreases = nnz(S.MeanDelta<0);
S.NumPointQuantileDecreases = nnz(S.QuantileDelta<0);
end
