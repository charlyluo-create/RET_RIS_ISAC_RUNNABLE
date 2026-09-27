function [Hd_user,Hr_user,Hd_tar,Hr_tar,G] = generate_channel( ...
    M,N,K,RIS2user_az,BS2tar_az,RIS2tar_az, ...
    loss_Bu,loss_Ru,loss_BR,loss_Bt,loss_Rt,beta_Ru, ...
    BS2RIS_az,RIS2BS_az)
%GENERATE_CHANNEL Geometry-consistent RET-RIS ISAC channels.
% Channel convention:
%   Hd_user(:,k), Hr_user(:,k), Hd_tar(:,t), Hr_tar(:,t) are column
%   channel vectors. For RIS phase phi, the effective communication column
%   channel is
%       h_k = h_d,k + G' * diag(conj(phi)) * h_r,k,
%   so h_k' * w is the received complex amplitude.
%
% The current implementation uses half-wavelength ULAs for the horizontal
% BS and RIS responses. BS-user direct links are Rayleigh; RIS-user links
% are Rician; BS-RIS and target links are LoS.

if nargin < 13 || isempty(BS2RIS_az); BS2RIS_az = 0; end
if nargin < 14 || isempty(RIS2BS_az); RIS2BS_az = 0; end

RIS2user_az = reshape(RIS2user_az,1,[]);
BS2tar_az = reshape(BS2tar_az,1,[]);
RIS2tar_az = reshape(RIS2tar_az,1,[]);
loss_Bu = reshape(loss_Bu,1,[]);
loss_Ru = reshape(loss_Ru,1,[]);
loss_Bt = reshape(loss_Bt,1,[]);
loss_Rt = reshape(loss_Rt,1,[]);

if numel(RIS2user_az) ~= K
    error('RIS2user_az must contain K=%d angles.',K);
end
T = numel(BS2tar_az);

%% BS -> user: NLoS Rayleigh
Hd_user = (randn(M,K)+1i*randn(M,K))/sqrt(2);
Hd_user = Hd_user .* sqrt(loss_Bu);

%% RIS -> user: Rician with independent LoS phase per user
nR = (0:N-1).';
a_RU = exp(-1i*pi*nR*sind(RIS2user_az));
los_phase = exp(1i*2*pi*rand(1,K));
Hr_LOS = a_RU .* los_phase;
Hr_NLOS = (randn(N,K)+1i*randn(N,K))/sqrt(2);
Hr_user = sqrt(loss_Ru) .* ( ...
    sqrt(beta_Ru/(1+beta_Ru))*Hr_LOS + ...
    sqrt(1/(1+beta_Ru))*Hr_NLOS);

%% BS -> RIS: rank-one geometric LoS MIMO channel
nB = (0:M-1).';
a_BS_BR = exp(-1i*pi*nB*sind(BS2RIS_az));
a_RIS_BR = exp(-1i*pi*nR*sind(RIS2BS_az));
G = sqrt(loss_BR) * a_RIS_BR * a_BS_BR';

%% BS -> target and RIS -> target: LoS steering vectors
Hd_tar = exp(-1i*pi*nB*sind(BS2tar_az)) .* sqrt(loss_Bt);
Hr_tar = exp(-1i*pi*nR*sind(RIS2tar_az)) .* sqrt(loss_Rt);

if size(Hd_tar,2) ~= T || size(Hr_tar,2) ~= T
    error('Target channel dimension mismatch.');
end
end
