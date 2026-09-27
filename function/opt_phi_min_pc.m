function phi = opt_phi_min_pc(N,K,M,Hd_user,Hr_user,G,W,A,rho,phi0)
%OPT_PHI_MIN_PC RCG minimization of the communication consistency penalty.
% Uses the same channel convention as build_effective_user_channel.m.
if nargin < 10 || isempty(phi0); phi0=ones(N,1); end
if size(W,1)~=M || size(W,2)~=K; error('W must be M-by-K.'); end
options=struct('maxiter',80,'tolgradnorm',1e-6,'minstepsize',1e-10, ...
    'restartperiod',12,'verbosity',0);
[phi,~,~]=unit_modulus_rcg(@costgrad,phi0,options);

    function [f,g]=costgrad(x)
        f=0; g=zeros(N,1);
        for kk=1:K
            hr=Hr_user(:,kk);
            for jj=1:K
                gw=G*W(:,jj);
                c=conj(hr).*gw;
                d=Hd_user(:,kk)'*W(:,jj)-A(kk,jj);
                s=d+c.'*x;
                f=f+rho*abs(s)^2;
                g=g+rho*2*conj(c)*s;
            end
        end
        f=real(f);
    end
end
