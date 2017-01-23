function definput=arg_adaptloop_(definput)
% function definput=arg_adaptloop_(definput)
%
% 1. Description:
%       Get the default values for the adaptation loop processing.
%       
%       'adt_dau1996', 'adt_dorp', 'adt_jepsen' added by AO
% 
% Author        : Peter L. Soendergaard and Piotr Majdak
% Downloaded on : 18/03/2014
% Modified by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Last update on: 22/04/2015 
% Last use on   : 01/07/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

definput.keyvals.limit  = 10;
definput.keyvals.minlvl =  0; 
definput.keyvals.tau    = [0.005 0.050 0.129 0.253 0.500]; % default
definput.keyvals.intnoise_dB = 0;

definput.groups.adt_dau       = {'tau',[0.005 0.050 0.129 0.253 0.500]}; % 0.005  0.050  0.129   0.253  0.500
definput.groups.adt_dau1996   = {'tau',[0.005 0.050 0.129 0.253 0.500],'limit',0}; % 0.005  0.050  0.129   0.253  0.500
definput.groups.adt_breebaart = {'tau',linspace(0.005,0.5,5)};                     % 0.005  0.1288 0.2525  0.3762 0.5000
definput.groups.adt_dorp      = {'tau',linspace(0.005,0.5,5),'limit',0};           % 0.005  0.1288 0.2525  0.3762 0.5000
definput.groups.adt_puschel   = {'tau',linspace(0.005,0.5,5),'limit',0};           % 0.005  0.1288 0.2525  0.3762 0.5000
definput.groups.adt_jepsen    = {'tau',[0.005 0.050 0.129 0.253 0.500],'limit',10,'intnoise_dB',9,'minlvl',-36}; % 0.005  0.050  0.129   0.253  0.500

%   Url: http://amtoolbox.sourceforge.net/doc/comp/arg_adaptloop.php

% Copyright (C) 2009-2015 Peter L. Soendergaard and Piotr Majdak.
% This file is part of AMToolbox version 0.9.5-0.9.7

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
