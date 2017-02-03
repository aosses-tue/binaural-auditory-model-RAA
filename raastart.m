function raastart
% function raastart
%
% 1. Description: 
%       - Initialises AMT Toolbox
%
% 2. Stand-alone example:
%       raastart;
%
% 4. Additional info:
%   Tested cross-platform: Yes
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2017
% Created on    : 20/01/2017
% Last edited on: 03/02/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dirmain  = [cd delim];
bSuccess = isdir([dirmain 'tb_AM_AddOns' delim]);

if bSuccess == 0
    error('Please set the current MATLAB directory to the location of the RAA model and re-run this script');
end

dir = [dirmain 'tb_AM_AddOns' delim]; addpath(dir);
subdirs = il_get_AM_AddOns_paths(dir);
for i = 1:length(subdirs)
  addpath(subdirs{i})
end
dir = [dirmain 'tb_AM' delim]; addpath(dir);
amtstart;

dir = [dirmain 'Utility'      delim]; addpath(dir);

disp(['EOF: ' mfilename])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function subdirs = il_get_AM_AddOns_paths(dir)

subdirs{1} = [dir 'arg'         delim];
subdirs{2} = [dir 'binaural'    delim];
subdirs{3} = [dir 'demos'       delim];
subdirs{4} = [dir 'filters'     delim];
subdirs{5} = [dir 'modelstages' delim];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function delimiter = delim

if isunix
    delimiter = '/';
else
    delimiter = '\';
end
