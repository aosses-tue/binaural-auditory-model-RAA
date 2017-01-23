function [outsig_dB, dBFS] = Do_SLM(insig,fs,weight_freq,weight_time,dBFS)
% function [outsig_dB, dBFS] = Do_SLM(insig,fs,weight_freq,weight_time,dBFS)
%
% 1. Description:
%       outsig is the weighted pressure in [Pa]
% 
% 2. Stand-alone example:
%   fs = 44100;
%   dur = 1; % s
%   insig = Create_sin(1000,dur,fs);
%   Do_SLM(insig,fs,'A','f');
% 
%   [insig, fs] = Wavread('D:\Databases\dir03-Speech\dutch\LISTman\jwz551.wav');
%   Do_SLM(insig,fs,'Z','f');
% 
% 3. Additional info:
%       Tested cross-platform: No
%       See also: Get_Leq
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 12/07/2016
% Last update on: 12/07/2016 
% Last use on   : 19/07/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    weight_freq = 'A';
end

if nargin < 4
    weight_time = 'f';
end

if nargin < 5
    dBFS = 100;
end

[b,a] = Gen_weighting_filters(fs,weight_freq);

% same processing, but without AMT, as done in PsySound:
dBoffset = 0.93; % determined empirically on 13/07/2016 to obtain the same values
                 % with this implementation and the one in PsySound.
calCoeff = 10.^((dBFS+dBoffset-94)/20);
insig = calCoeff*insig; % same as: insig = setdbspl(insig,rmsdb(insig)+dBFS,'dboffset',94);

outsig = filter(b,a,insig);
outsig = il_integrator(outsig,fs,weight_time);

outsig_dB = 20*log10(abs(outsig)/2e-5);

if nargout == 0
    Leq = Get_Leq(outsig_dB);
    figure;
    plot( outsig_dB ); grid on;
    title(sprintf('Level Leq=%.1f [dB]',Leq));
    disp('')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EOF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outsig = il_integrator(insig,fs, weightingType)
% function outsig = il_integrator(insig,fs, weightingType)
% 
% INTEGRATOR 
% 
% Generates and applies the integration filter to the input data
%
% For the theory behind this code, please refer to the document 'psy-sound.pdf'
%
% Author : Matt R. Flax <flatmax>
%          Jan. 2007 for the psysound project.
%
% Revised : Farhan Rizwi
%           Mostly cleanup of redundant code
%           July 2007 for the psysound project.
%
% input :
%         dataIn     - data vector
%         fastOrSlow - RC time constant is 'f' fast (125 ms) or 's'
%                      slow (1 s)
%         Fs         - sample rate of the data
%     
% Time constant for the leaky integrator is tau. This is basically
% a low-pass filter with a bandwidth of 1/tau.
%
% This yields the following transfer function :
%                      1
%        H(s) =  ---------------
%                tau s  +   1
%
% This can also be used with an rms integrator - 
% make fastorslow = 'r0.2' for an rms integrator with a window size of 0.2
% seconds. For a leaky integrator or 0.2 seconds 'l0.2' will work. 

% Filter coeffecients

switch weightingType
    case 'f' % fast leak - time constant = 125 ms
        tau = 125e-3;
    case 's' % slow leak - time constant = 1 s
        tau = 1;
    case 'i'
        tau = 35e-3; % impulse
    case 'p'
        tau = 50e-6;	
    case 'l'
        tau = str2num(fastOrSlow(2:end));	
    case 'r'
        tau = str2num(fastOrSlow(2:end));	
    otherwise
        error(['integrator: unknown leak case ' char(fastOrSlow)]);
end

% State vector
Z = [];

% Exponential term
E = exp(-1/(tau*fs));

% Filter numerator - with gain adjustment
b = 1 - E;

% Filter denominator
a = [1 -E];

% State vector
Z = [];

% Create run function handle
% Use filter to perform the integration
[outsig, Z] = filter(b, a, abs(insig), Z, 1);
