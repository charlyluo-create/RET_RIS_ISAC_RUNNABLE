function [Wr,p_radar] = allocate_radar_power(Hu,C_radar,P_total,p_comm)
%ALLOCATE_RADAR_POWER Allocate residual power in Null(Hu').
M=size(Hu,1); Wr=zeros(M,1);
p_radar=max(real(P_total-p_comm),0);
if p_radar<=1e-12 || ~isfinite(p_radar); p_radar=0; return; end
Z=null(Hu');
if isempty(Z); p_radar=0; return; end
B=Z'*C_radar*Z; B=(B+B')/2;
[V,D]=eig(B); [~,idx]=max(real(diag(D)));
v=Z*V(:,idx);
if norm(v)<=1e-12; v=Z(:,1); end
Wr=sqrt(p_radar)*v/norm(v);
end
