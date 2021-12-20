function [pRev,Level,outputs] = osses2017(filename,dBFS)
% function [pRev,Level,outputs] = osses2017(filename,dBFS)
%
% 1. Description:
%       This script implements the Room Acoustic Analyser (RAA) model to
%       compute reverberance estimates pRev from an audio recording (filename)
% 
%    Input parameters:
%       filename - file name of the sound to be processed. Make sure that 
%                  the file is 'visible' for MATLAB. The sound should ideally
%                  be a stereo Wav file, but in case of a mono signal a diotic
%                  condition is assumed (left-ear signal is copied to the right-
%                  ear one)
%       dBFS     - Calibration level to convert dBFS to dB SPL. By default an
%                  amplitude of 1 (in the waveform) is assumed to represent 
%                  a level of 100 dB SPL (default in AMT standard)
% 
%    Output parameters:
%       pRev - reverberance estimates in Model Units
%           pRev(1) = median pRev
%           pRev(2) = minimum pRev
%           pRev(3) = maximum pRev
%           pRev(4) = N analysis frames
%       Level - array with overall levels dB:
%           Level(1) is Leq in dB(A)
%           Level(2) is Lmax in dB(A)
%           Level(3) is Leq in dB linear or dB(Z)
%           Level(4) is Lmax in dB linear or dB(Z)
%       outputs - struct containing additional information
%           outputs.pClar - provides the estimate of perceptual clarity
% 
% 2. Stand-alone example:
%       % Process file: b06-Cello-room-5-B.wav
%       dBFS = 119; % amplitude 1 = 119 dB SPL (information from recording setup)
%       dir_where = raabasepath; % Example path in Windows:  dir_where = 'D:\MATLAB_RAA\tb_AM_AddOns\';
%                                % Example path in Unix sys: dir_where = '/home/alejandro/Documenten/MATLAB_RAA/tb_AM_AddOns/';
%       filename = [dir_where 'auxdata' filesep 'osses2017' filesep 'b06-Cello-room-5-B.wav'];
%       [pRev1,Level1] = demo_raa(filename,dBFS);
%       % Expected result:
%       %   pRev1  = 18.7360   16.7706   19.8998    6.0000
%       %   Level1 = 66.2660   72.6133   73.8960   82.0570
% 
% 3. Additional info:
%       See reference Osses2017
%       Tested cross-platform: Yes
% 
% Old name in RAA v1.0: demo_raa.m
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2017
% Created on    : 21/01/2017
% Last update on: 03/02/2017 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    dBFS = dbspl(1); % 100 dB SPL attributed to an amplitude of 1, AMT default
end

Level = zeros(1,4); % Initialisation

[insig,fs] = audioread(filename);

subfs = 11025; % Hz

insig = insig(1:round(10*fs),:);

insig = gaindb(1,dBFS-100)*insig; % gaindb replaces 'my' From_dB

dB = Do_SLM(insig,fs,'A','f',100);
if size(dB,2) == 2
    % If signal is stereo, the SPL are averaged into one channel
    dB = 10*log10(0.5*abs(10.^(dB(:,1)/10)+10.^(dB(:,2)/10)));
end
    
Level(1) = Get_Leq(dB); % LAeq
Level(2) = max(dB); % LAmax

dB = Do_SLM(insig,fs,'Z','f',100);
if size(dB,2) == 2
    % If signal is stereo, the SPL are averaged into one channel
    dB = 10*log10(0.5*abs(10.^(dB(:,1)/10)+10.^(dB(:,2)/10)));
end
Level(3) = Get_Leq(dB); % LZeq
Level(4) = max(dB); % LZmax

[outsig, fc, par] = dorp2011(insig, fs, 'subfs',subfs);

pRev(1) = median(par.pRev_frame);
pRev(2) = min(par.pRev_frame);
pRev(3) = max(par.pRev_frame);
pRev(4) = length(par.IdxFrames); % N frames

pClar(1) = median(par.pClar_frame);
pClar(2) = min(par.pClar_frame);
pClar(3) = max(par.pClar_frame);
pClar(4) = length(par.IdxFrames); % N frames

par = rmfield(par,'t_psi');
if nargout >= 3
    outputs.pClar = pClar;
    outputs.par = par;
end
