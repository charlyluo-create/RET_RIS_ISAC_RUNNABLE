function C_radar = build_radar_matrix(use_ris,phi,Hd_tar,Hr_tar,G,sigma_r2,omega)
%BUILD_RADAR_MATRIX Weighted monostatic four-path sensing matrix.
% One-way target column channel:
%   h_t = h_Bt + G' diag(conj(phi)) h_Rt.
% Under channel reciprocity, the round-trip matrix is
%   H_t = conj(h_t) h_t^T,
% whose expansion contains direct-direct, direct-RIS, RIS-direct and
% RIS-RIS paths. The reported metric is trace(W'*C_radar*W).
[M,T]=size(Hd_tar);
omega=omega(:); omega=omega/max(sum(omega),eps);
C_radar=zeros(M,M);
for t=1:T
    if use_ris
        h=Hd_tar(:,t)+G'*(conj(phi(:)).*Hr_tar(:,t));
    else
        h=Hd_tar(:,t);
    end
    Ht=conj(h)*h.';
    C_radar=C_radar+omega(t)*(Ht'*Ht)/max(sigma_r2,eps);
end
C_radar=(C_radar+C_radar')/2;
end
