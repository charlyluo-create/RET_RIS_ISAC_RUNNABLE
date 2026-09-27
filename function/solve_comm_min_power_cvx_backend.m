function [Wc,p_comm,feasible,status_out,meta] = solve_comm_min_power_cvx_backend(Hu,gamma_vec,sigma_c)
%SOLVE_COMM_MIN_POWER_CVX_BACKEND CVX/SOCP reference solver for QoS beamforming.
%
% This function is intentionally isolated from telemetry/nested callbacks.
% CVX creates declared variables dynamically in the caller workspace; keeping
% this workspace non-static avoids MATLAB R2025b errors such as:
%   "Invalid syntax for calling function 'X' on the path ..."
%
% Problem:
%   min_W ||W||_F^2
%   s.t. SINR_k >= gamma_k,  k=1,...,K.
%
% The standard phase rotation converts each SINR constraint to an SOC form.

[M,K]=size(Hu);
Wc=zeros(M,K); p_comm=Inf; feasible=false; status_out='CVX not started';
meta=struct('min_sinr_ratio',NaN,'cvx_status','Not started');

gamma_vec=gamma_vec(:);
if isempty(Hu) || numel(gamma_vec)~=K || sigma_c<0 || ...
        any(gamma_vec<=0) || any(~isfinite(gamma_vec)) || ~all(isfinite(Hu(:)))
    status_out='CVX invalid input';
    return;
end

% Numerical scaling leaves every SINR unchanged while improving conditioning.
scale=max(norm(Hu,'fro'),1e-12);
Hs=Hu/scale;
ns=sigma_c/scale;

try; cvx_clear; catch; end

% Explicit initialization prevents MATLAB from resolving the symbol as a
% path function before CVX replaces it with a cvx variable.
Wcvx=[]; %#ok<NASGU>

try
    cvx_begin quiet
        variable Wcvx(M,K) complex
        minimize( square_pos(norm(Wcvx,'fro')) )
        subject to
            for k=1:K
                hk=Hs(:,k)';
                imag(hk*Wcvx(:,k)) == 0;
                real(hk*Wcvx(:,k)) >= 0;
                norm([hk*Wcvx, ns],2) <= ...
                    sqrt(1+1/gamma_vec(k))*real(hk*Wcvx(:,k));
            end
    cvx_end

    meta.cvx_status=cvx_status;
    status_out=cvx_status;

    if contains(cvx_status,'Solved') && all(isfinite(Wcvx(:)))
        sinr=zeros(K,1);
        for k=1:K
            y=Hu(:,k)'*Wcvx;
            denom=sum(abs(y).^2)-abs(y(k))^2+sigma_c^2;
            sinr(k)=abs(y(k))^2/max(real(denom),realmin);
        end
        meta.min_sinr_ratio=min(sinr./gamma_vec);

        if all(sinr>=gamma_vec*(1-5e-4))
            Wc=Wcvx;
            p_comm=real(norm(Wcvx,'fro')^2);
            feasible=isfinite(p_comm) && p_comm>=0;
        else
            status_out=[cvx_status ' (SINR validation failed)'];
        end
    end
catch ME
    status_out=['CVX error: ' ME.message];
    meta.cvx_status=status_out;
    try; cvx_clear; catch; end
end
end
