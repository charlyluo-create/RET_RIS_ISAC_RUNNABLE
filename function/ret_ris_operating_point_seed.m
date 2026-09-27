function seed = ret_ris_operating_point_seed(base,gamma,P_total)
%RET_RIS_OPERATING_POINT_SEED Deterministic seed keyed to (gamma,P_total).
% This makes an overlapping operating point identical in Exp. 2 and Exp. 3.
gamma_dB=10*log10(max(gamma,realmin));
P_dBm=10*log10(max(P_total,realmin)/1e-3);
gkey=round((gamma_dB+100)*1000);
pkey=round((P_dBm+200)*1000);
modulus=2^32-1;
seed=mod(double(base)+73856093*double(gkey)+19349663*double(pkey),modulus);
seed=floor(seed); if seed<0; seed=seed+modulus; end
end
