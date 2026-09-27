function scen = apply_relative_motion_phase(scen,scen_ref,lambda)
%APPLY_RELATIVE_MOTION_PHASE Apply paired path-length phase evolution.
%
% The same normalized fading realization is used at the reference and moved
% snapshots. This helper adds only the deterministic carrier-phase rotation
% caused by the path-length difference, so the mobility experiment does not
% confuse channel aging with an independent fading redraw.
%
% At the reference snapshot all phase rotations equal one.

if nargin < 3 || ~isscalar(lambda) || ~isfinite(lambda) || lambda <= 0
    error('A positive carrier wavelength is required.');
end
required={'d_Bu','d_Ru','d_Bt','d_Rt','Hd_user','Hr_user','Hd_tar','Hr_tar'};
for ii=1:numel(required)
    if ~isfield(scen,required{ii}) || ~isfield(scen_ref,required{ii})
        error('Scenario field %s is required for paired phase evolution.',required{ii});
    end
end

phase_Bu=exp(-1i*2*pi*(scen.d_Bu-scen_ref.d_Bu)/lambda);
phase_Ru=exp(-1i*2*pi*(scen.d_Ru-scen_ref.d_Ru)/lambda);
phase_Bt=exp(-1i*2*pi*(scen.d_Bt-scen_ref.d_Bt)/lambda);
phase_Rt=exp(-1i*2*pi*(scen.d_Rt-scen_ref.d_Rt)/lambda);

scen.Hd_user=scen.Hd_user.*reshape(phase_Bu,1,[]);
scen.Hr_user=scen.Hr_user.*reshape(phase_Ru,1,[]);
scen.Hd_tar=scen.Hd_tar.*reshape(phase_Bt,1,[]);
scen.Hr_tar=scen.Hr_tar.*reshape(phase_Rt,1,[]);
end
