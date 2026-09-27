function policy = get_figure_policy(run_mode)
%GET_FIGURE_POLICY Centralized figure-output policy for all experiments.
% Environment variable RET_RIS_FIGURE_PROFILE accepts:
%   paper_master : save no per-experiment figures; the master figure builder
%                  creates the six composite manuscript figures after all data
%                  files have been produced.
%   compact      : save only one primary/composite figure per experiment.
%   diagnostic   : save the primary figure and every supplementary diagnostic.
%   none         : save no figures (useful for FAST validation).
%
% When the environment variable is absent, PAPER runs default to compact and
% FAST runs default to none. Computed MAT/CSV results are unaffected.

if nargin < 1 || isempty(run_mode)
    run_mode = 'fast';
end
profile = lower(strtrim(getenv('RET_RIS_FIGURE_PROFILE')));
if isempty(profile)
    if strcmpi(run_mode,'paper')
        profile = 'compact';
    else
        profile = 'none';
    end
end
valid = {'paper_master','compact','diagnostic','none'};
if ~ismember(profile,valid)
    warning('Unknown RET_RIS_FIGURE_PROFILE="%s". Using compact.',profile);
    profile = 'compact';
end
policy = struct();
policy.profile = profile;
policy.save_main = ismember(profile,{'compact','diagnostic'});
policy.save_supplementary = strcmp(profile,'diagnostic');
policy.save_any = policy.save_main || policy.save_supplementary;
policy.build_master = strcmp(profile,'paper_master');
policy.silent = ismember(profile,{'paper_master','none'});
end
