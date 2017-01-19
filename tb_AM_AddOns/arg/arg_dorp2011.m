function definput = arg_dorp2011(definput)
% function definput = arg_dorp2011(definput)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 11/07/2016
% Last update on: 29/07/2016 
% Last use on   : 18/11/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% General options of the model:
definput.keyvals.order  = 1; % order outer- middle-ear filter
% Set option outermiddleear to 'outermiddleear' and therefore: do_outermiddleear to 1
definput.flags.outermiddleear={'outermiddleear','nooutermiddleear'}; 
definput.flags.binaural_proc = {'binaural','no_binaural'};

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
