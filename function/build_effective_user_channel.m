function Hu = build_effective_user_channel(phi,use_ris,Hd_user,Hr_user,G)
%BUILD_EFFECTIVE_USER_CHANNEL Return M-by-K effective column channels.
% h_k = h_d,k + G' diag(conj(phi)) h_r,k.
if ~use_ris
    Hu = Hd_user;
    return;
end
N = size(G,1);
phi = phi(:);
if numel(phi) ~= N
    error('RIS phase length (%d) must equal N=%d.',numel(phi),N);
end
bad = abs(phi)<1e-12; phi(bad)=1; phi=phi./abs(phi);
Hu = Hd_user + G' * (conj(phi).*Hr_user);
end
