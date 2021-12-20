function bp = raabasepath
% function bp = raabasepath
%
% 1. Description:
%       Base path of the RAA binaural model relative to the AMT_AddOns 
%       directory.
% 
%   Usage: bp = amtbasepath;
%
%   `amtbasepath` returns the top level directory in which the AMT
%   files are installed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
bp = fileparts(fileparts(mfilename('fullpath')));
bp = [bp filesep]; % Adding ('/' or '\')