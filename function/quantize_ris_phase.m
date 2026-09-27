function phi_q = quantize_ris_phase(phi,bits)
%QUANTIZE_RIS_PHASE Nearest-neighbour uniform phase quantization.
% BITS=Inf returns the normalized continuous phase vector.

phi=phi(:);
if isempty(phi)
    phi_q=phi;
    return;
end
bad=abs(phi)<1e-12;
phi(bad)=1;
phi=phi./abs(phi);
if isinf(bits)
    phi_q=phi;
    return;
end
if ~isscalar(bits) || bits<1 || abs(bits-round(bits))>1e-12
    error('bits must be a positive integer or Inf.');
end
L=2^round(bits);
ang=mod(angle(phi),2*pi);
idx=mod(round(L*ang/(2*pi)),L);
phi_q=exp(1i*2*pi*idx/L);
end
