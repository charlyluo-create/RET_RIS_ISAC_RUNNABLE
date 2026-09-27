function scen = generate_isac_scenario(cfg,N,K,seed,x_u,y_u,common_direct)
%GENERATE_ISAC_SCENARIO Generate one unified roadside RET-RIS ISAC sample.
%
% Inputs
%   cfg           : unified configuration structure
%   N, K          : RIS element count and communication-user count
%   seed          : deterministic random seed
%   x_u, y_u      : optional 1-by-K user coordinates. Empty = random road users
%   common_direct : optional structure with Hd_user and Hd_tar. It is useful
%                   for paired scalability comparisons across different N.
%
% Output
%   scen contains channels, RET gain tables, geometry and useful tilt indices.

% Use one deterministic random stream. When coordinates are not supplied,
% x and y are drawn sequentially exactly as in Experiments 1--3.
rng(seed,'twister');
if nargin < 5 || isempty(x_u)
    x_u = cfg.road_x_min + (cfg.road_x_max-cfg.road_x_min)*rand(1,K);
end
if nargin < 6 || isempty(y_u)
    y_u = -cfg.road_y_half + 2*cfg.road_y_half*rand(1,K);
end
if nargin < 7 || isempty(common_direct)
    common_direct = struct();
end
% Optional target-position override is used by the mobility experiment while
% preserving the original calling convention for Experiments 1--5 and 7.
target_x = cfg.target_x;
target_y = cfg.target_y;
if isfield(common_direct,'target_x') && ~isempty(common_direct.target_x)
    target_x = reshape(common_direct.target_x,1,[]);
end
if isfield(common_direct,'target_y') && ~isempty(common_direct.target_y)
    target_y = reshape(common_direct.target_y,1,[]);
end
if numel(target_x) ~= cfg.T || numel(target_y) ~= cfg.T
    error('Target-position override must contain T=%d entries.',cfg.T);
end
x_u = reshape(x_u,1,[]);
y_u = reshape(y_u,1,[]);
if numel(x_u) ~= K || numel(y_u) ~= K
    error('x_u and y_u must each contain K=%d entries.',K);
end

M=cfg.M; T=cfg.T;
theta_range=cfg.theta_tilt_range;
numTheta=numel(theta_range);

% Horizontal geometry controls azimuth and electrical downtilt.
d_Bu_h=hypot(x_u-cfg.bs_x,y_u-cfg.bs_y);
d_Ru_h=hypot(x_u-cfg.ris_x,y_u-cfg.ris_y);
RIS2user_az=atan2d(y_u-cfg.ris_y,x_u-cfg.ris_x);

% Three-dimensional propagation distance controls path loss.
d_Bu=sqrt(d_Bu_h.^2+(cfg.h_BS-cfg.h_user).^2);
d_Ru=sqrt(d_Ru_h.^2+(cfg.h_RIS-cfg.h_user).^2);

loss_BR=cfg.C0*(cfg.d_BR/cfg.D0)^(-cfg.alpha_BR);
loss_Bu=10^(-cfg.direct_extra_loss_dB/10)* ...
    cfg.C0*(d_Bu/cfg.D0).^(-cfg.alpha_Bu);
loss_Ru=cfg.C0*(d_Ru/cfg.D0).^(-cfg.alpha_Ru);
d_Bt_h=hypot(target_x-cfg.bs_x,target_y-cfg.bs_y);
d_Rt_h=hypot(target_x-cfg.ris_x,target_y-cfg.ris_y);
d_Bt=sqrt(d_Bt_h.^2+(cfg.h_BS-cfg.h_target).^2);
d_Rt=sqrt(d_Rt_h.^2+(cfg.h_RIS-cfg.h_target).^2);
target_az=atan2d(target_y-cfg.bs_y,target_x-cfg.bs_x);
RIS2target_az=atan2d(target_y-cfg.ris_y,target_x-cfg.ris_x);
loss_Bt=cfg.C0*(d_Bt/cfg.D0).^(-cfg.alpha_Bt);
loss_Rt=cfg.C0*(d_Rt/cfg.D0).^(-cfg.alpha_Rt);

[Hd_user,Hr_user,Hd_tar,Hr_tar,G]=generate_channel( ...
    M,N,K,RIS2user_az,target_az,RIS2target_az, ...
    loss_Bu,loss_Ru,loss_BR,loss_Bt,loss_Rt,cfg.beta_Ru, ...
    cfg.bs2ris_az_deg,cfg.ris2bs_az_deg);

if isfield(common_direct,'Hd_user') && ~isempty(common_direct.Hd_user)
    if size(common_direct.Hd_user,1)~=M || size(common_direct.Hd_user,2)<K
        error('common_direct.Hd_user must be M-by-at-least-K.');
    end
    Hd_user=common_direct.Hd_user(:,1:K);
end
if isfield(common_direct,'Hd_tar') && ~isempty(common_direct.Hd_tar)
    if ~isequal(size(common_direct.Hd_tar),[M T])
        error('common_direct.Hd_tar must be M-by-T.');
    end
    Hd_tar=common_direct.Hd_tar;
end

% RET gains for direct users, BS-RIS and direct target illumination.
elev_users=atan2d(cfg.h_BS-cfg.h_user,d_Bu_h(:));
elev_RIS=atan2d(cfg.h_BS-cfg.h_RIS,cfg.d_BR_horizontal);
elev_targets=atan2d(cfg.h_BS-cfg.h_target,d_Bt_h(:));

gain_u=zeros(numTheta,K);
gain_ris=zeros(numTheta,1);
gain_t=zeros(numTheta,T);
for ti=1:numTheta
    for k=1:K
        gain_u(ti,k)=ret_vertical_gain(elev_users(k),theta_range(ti), ...
            cfg.theta_3dB,cfg.SLA_V,cfg.N_RET,cfg.d_RET_lambda);
    end
    gain_ris(ti)=ret_vertical_gain(elev_RIS,theta_range(ti), ...
        cfg.theta_3dB,cfg.SLA_V,cfg.N_RET,cfg.d_RET_lambda);
    for t=1:T
        gain_t(ti,t)=ret_vertical_gain(elev_targets(t),theta_range(ti), ...
            cfg.theta_3dB,cfg.SLA_V,cfg.N_RET,cfg.d_RET_lambda);
    end
end

[~,ref_idx]=min(abs(theta_range-cfg.theta_ref));
target_elev_mean=sum(cfg.target_weights(:).*elev_targets(:));
[~,target_idx]=min(abs(theta_range-target_elev_mean));

scen=struct();
scen.M=M; scen.N=N; scen.K=K; scen.T=T;
scen.x_u=x_u; scen.y_u=y_u;
scen.target_x=target_x; scen.target_y=target_y;
scen.d_Bu=d_Bu; scen.d_Ru=d_Ru; scen.d_Bt=d_Bt; scen.d_Rt=d_Rt;
scen.RIS2user_az=RIS2user_az;
scen.target_az=target_az; scen.RIS2target_az=RIS2target_az;
scen.Hd_user=Hd_user; scen.Hr_user=Hr_user;
scen.Hd_tar=Hd_tar; scen.Hr_tar=Hr_tar; scen.G=G;
scen.gain_u=gain_u; scen.gain_ris=gain_ris; scen.gain_t=gain_t;
scen.theta_range=theta_range;
scen.ref_idx=ref_idx; scen.target_idx=target_idx;
scen.elev_users=elev_users; scen.elev_targets=elev_targets;
end
