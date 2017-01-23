function definput = arg_dorp2011(definput)
% function definput = arg_dorp2011(definput)
%
% 1. Description:
%       Loads arguments to be used by the dorp2011preproc function.
% 
% 2. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2017
% Created on    : 11/07/2016
% Last update on: 23/01/2017 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% General options of the model:
definput.keyvals.order  = 1; % order outer- middle-ear filter
definput.flags.outermiddleear={'outermiddleear','nooutermiddleear'}; 
definput.flags.binaural_proc = {'no_binaural','binaural'}; % binaural processor not implemented yet

%% For binaural processor:
definput.keyvals.muASW = 2.00e-2;
definput.keyvals.nuASW = 5.63e2;
definput.keyvals.muLEV = 2.76e-2;
definput.keyvals.nuLEV = 6.80e2;

%% For central processor:
definput.keyvals.Psimin = 0.34; % 'Optimisation' by running m20160714_fourth_report_MBBM_testing_vanDorp
                                % criterion used: minimise error (mean absolute difference).
% Psimin_dip determined from the ratio abs(-7.49e-3/1.33e-3), which are the
% Psimin and Psimin_dip as suggested by van Dorp:
definput.keyvals.Psimin_dip = -definput.keyvals.Psimin/5.63; 
definput.keyvals.Tmin   = 63.1e-3; % 63.1 [ms]
definput.keyvals.framelen = 5; % s
definput.keyvals.hopsize  = 1; % s

definput.flags.multiple_one_second={'multiple_one_second','not_multiple'};
definput.flags.exclude_silence={'exclude','not_exclude'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
