function cfg = refresh_ret_ris_geometry(cfg)
%REFRESH_RET_RIS_GEOMETRY Recompute all geometry-dependent config fields.
% Call this helper after changing BS, RIS, road-centre, or target positions.

required={'bs_x','bs_y','h_BS','ris_x','ris_y','h_RIS', ...
    'h_user','h_target','target_x','target_y','road_x_min','road_x_max', ...
    'ret_angle_resolution_deg'};
for ii=1:numel(required)
    if ~isfield(cfg,required{ii})
        error('Missing configuration field: %s',required{ii});
    end
end

cfg.d_BR_horizontal=hypot(cfg.ris_x-cfg.bs_x,cfg.ris_y-cfg.bs_y);
cfg.bs2ris_az_deg=atan2d(cfg.ris_y-cfg.bs_y,cfg.ris_x-cfg.bs_x);
cfg.ris2bs_az_deg=atan2d(cfg.bs_y-cfg.ris_y,cfg.bs_x-cfg.ris_x);

cfg.d_Bt_horizontal=hypot(cfg.target_x-cfg.bs_x,cfg.target_y-cfg.bs_y);
cfg.d_Rt_horizontal=hypot(cfg.target_x-cfg.ris_x,cfg.target_y-cfg.ris_y);
cfg.target_az_deg=atan2d(cfg.target_y-cfg.bs_y,cfg.target_x-cfg.bs_x);
cfg.RIS2tar_deg=atan2d(cfg.target_y-cfg.ris_y,cfg.target_x-cfg.ris_x);

cfg.d_BR=sqrt(cfg.d_BR_horizontal.^2+(cfg.h_BS-cfg.h_RIS).^2);
cfg.d_Bt_vec=sqrt(cfg.d_Bt_horizontal.^2+(cfg.h_BS-cfg.h_target).^2);
cfg.d_Rt_vec=sqrt(cfg.d_Rt_horizontal.^2+(cfg.h_RIS-cfg.h_target).^2);

cfg.fixed_ret_reference_x=0.5*(cfg.road_x_min+cfg.road_x_max);
cfg.fixed_ret_reference_y=0;
cfg.fixed_ret_raw_deg=atan2d(cfg.h_BS-cfg.h_user, ...
    hypot(cfg.fixed_ret_reference_x-cfg.bs_x, ...
          cfg.fixed_ret_reference_y-cfg.bs_y));
res=max(cfg.ret_angle_resolution_deg,eps);
cfg.theta_ref=res*round(cfg.fixed_ret_raw_deg/res);
end
