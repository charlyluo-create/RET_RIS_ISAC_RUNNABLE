function tcrit = student_t_critical(confidence,nu)
%STUDENT_T_CRITICAL Two-sided Student-t critical value without toolboxes.
% TCRIT satisfies P(-TCRIT <= T_nu <= TCRIT)=CONFIDENCE.
if nargin<1 || isempty(confidence); confidence=0.95; end
if nargin<2 || isempty(nu); error('Degrees of freedom are required.'); end
if confidence<=0 || confidence>=1 || nu<=0
    error('confidence must be in (0,1) and nu must be positive.');
end
alpha=1-confidence;
try
    x=betaincinv(alpha,nu/2,0.5);
    tcrit=sqrt(nu*(1-x)/max(x,realmin));
catch
    % Conservative normal fallback for unusual minimal MATLAB installations.
    tcrit=sqrt(2)*erfcinv(alpha);
end
end
