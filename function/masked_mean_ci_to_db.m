function [mu_db,lo_db,hi_db,n,mu_linear] = masked_mean_ci_to_db(X,mask,reference)
%MASKED_MEAN_CI_TO_DB Mean/95% CI in the linear domain, displayed in dB.
% REFERENCE is 1 for dimensionless SNR/gain and 1e-3 for W-to-dBm.
if nargin<3 || isempty(reference); reference=1; end
[mu,lo,hi,n]=masked_mean_ci_linear(X,mask);
mu_linear=mu;
mu_db=nan(size(mu)); lo_db=nan(size(mu)); hi_db=nan(size(mu));
valid=isfinite(mu) & mu>0;
mu_db(valid)=10*log10(mu(valid)/reference);
valid=isfinite(lo);
lo_db(valid)=10*log10(max(lo(valid),realmin)/reference);
valid=isfinite(hi);
hi_db(valid)=10*log10(max(hi(valid),realmin)/reference);
end
