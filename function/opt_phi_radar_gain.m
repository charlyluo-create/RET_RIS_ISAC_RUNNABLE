function phi = opt_phi_radar_gain(N,K,M,T,Hd_tar,Hr_tar,Hd_user,Hr_user,G,W,A,rho,omega,phi0)
%OPT_PHI_RADAR_GAIN Riemannian conjugate-gradient phase update.
% Minimizes -weighted one-way target gain plus the communication consistency
% penalty. Outer scripts accept a candidate only when the exact four-path
% radar objective improves.
if nargin < 14 || isempty(phi0); phi0=ones(N,1); end
omega=omega(:); omega=omega/max(sum(omega),eps);
options=struct('maxiter',80,'tolgradnorm',1e-6,'minstepsize',1e-10, ...
    'restartperiod',12,'verbosity',0);
[phi,~,~]=unit_modulus_rcg(@costgrad,phi0,options);

    function [f,g]=costgrad(x)
        radar=0; grad_r=zeros(N,1);
        for tt=1:T
            hr=Hr_tar(:,tt);
            D=conj(hr).*G; % N-by-M, implicit row scaling
            s=Hd_tar(:,tt)' + x.'*D; % 1-by-M effective row channel
            radar=radar+omega(tt)*real(s*s');
            grad_r=grad_r+omega(tt)*2*conj(D)*s.';
        end

        penalty=0; grad_p=zeros(N,1);
        for kk=1:K
            hr=Hr_user(:,kk);
            for jj=1:size(W,2)
                gw=G*W(:,jj);
                c=conj(hr).*gw;
                if kk<=size(A,1) && jj<=size(A,2)
                    aij=A(kk,jj);
                else
                    aij=0;
                end
                d=Hd_user(:,kk)'*W(:,jj)-aij;
                z=d+c.'*x;
                penalty=penalty+abs(z)^2;
                grad_p=grad_p+2*conj(c)*z;
            end
        end
        f=real(-radar+rho*penalty);
        g=-grad_r+rho*grad_p;
    end
end
