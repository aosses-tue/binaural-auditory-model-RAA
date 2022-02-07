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

dirmain  = [cd filesep];
bSuccess = exist([dirmain 'tb_AM_AddOns' filesep],'dir');

if bSuccess == 0
    error('Please set the current MATLAB directory to the location of the RAA model and re-run this script');
end

dir = [dirmain 'tb_AM_AddOns' filesep]; addpath(dir);
subdirs = il_get_AM_AddOns_paths(dir);
for i = 1:length(subdirs)
  addpath(subdirs{i})
end
if exist('amt_start.m','file')
    amt_start;
    % dir = [dirmain 'tb_AM' filesep]; addpath(dir);
    % amt_start;
else
    fprintf('Please download the latest version of AMT toolbox (amtoolbox.org/) and run amt_start before starting to use the RAA model\n');
end

dir = [dirmain 'Utility' filesep]; addpath(dir);

disp(['EOF: ' mfilename])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function subdirs = il_get_AM_AddOns_paths(dir)

subdirs{1} = [dir 'arg'         filesep];
subdirs{2} = [dir 'binaural'    filesep];
subdirs{3} = [dir 'common'      filesep];
subdirs{4} = [dir 'experiments' filesep];
subdirs{5} = [dir 'models'      filesep];
subdirs{6} = [dir 'modelstages' filesep];
