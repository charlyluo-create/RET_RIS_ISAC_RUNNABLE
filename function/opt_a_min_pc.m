function A = opt_a_min_pc(K,M,A,Hu,W,rho,sigma_c,gamma)
%OPT_A_MIN_PC Auxiliary-variable update for the communication SINR penalty.
A_prev=A;
try; cvx_clear; catch; end
try
    cvx_begin quiet
        variable X(K,K) complex
        expression obj
        obj=0;
        for k=1:K
            for j=1:K
                obj=obj+square_abs(X(k,j)-Hu(:,k)'*W(:,j));
            end
        end
        minimize(rho*obj)
        subject to
            for k=1:K
                imag(X(k,k))==0;
                real(X(k,k))>=0;
                norm([X(k,:),sigma_c],2) <= sqrt(1+1/gamma)*real(X(k,k));
            end
    cvx_end
    if contains(cvx_status,'Solved') && all(isfinite(X(:)))
        A=X;
    else
        A=A_prev;
    end
catch
    A=A_prev;
    try; cvx_clear; catch; end
end
end
