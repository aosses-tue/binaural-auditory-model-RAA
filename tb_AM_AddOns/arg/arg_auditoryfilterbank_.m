function definput=arg_auditoryfilterbank_(definput)
% function definput=arg_auditoryfilterbank_(definput)
%
% 1. Description:
%       Default parameters for Gammatone filter bank.
% 
%       The config 'gtf_dorp2011' needs to be added to AMT 1.0 (added by AO)
% 
% Last used on: 27/06/2015
% Last edited on: 21/12/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

definput.keyvals.flow=80;
definput.keyvals.fhigh=8000;
definput.keyvals.basef=[];
definput.keyvals.bwmul=1;

definput.groups.gtf_dorp2011 = {'flow',168,'fhigh',1840};
