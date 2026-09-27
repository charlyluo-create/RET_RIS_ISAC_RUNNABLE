function gain = ret_vertical_gain(elev_deg,theta_tilt_deg,theta_3dB,SLA_V,N_RET,d_lambda)
%RET_VERTICAL_GAIN Sector-common electronic downtilt gain for the BS panel.
% One common theta_tilt_deg is shared by all active antenna ports. Angles are
% downward from the local horizontal. The unit-norm RET weights
% redirect spatial power without changing the digital beamforming power.
A_EV_dB = -min(12*(elev_deg/theta_3dB).^2,SLA_V);
G_elem = 10.^(A_EV_dB/10);
n = (0:N_RET-1).';
a = exp(-1i*2*pi*d_lambda*n*sind(elev_deg));
w = exp(-1i*2*pi*d_lambda*n*sind(theta_tilt_deg))/sqrt(N_RET);
G_array = abs(w'*a).^2;
gain = max(real(G_elem.*G_array),1e-12);
end
