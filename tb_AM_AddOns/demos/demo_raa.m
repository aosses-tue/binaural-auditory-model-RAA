function [pRev,Level,outputs] = demo_raa(filename,dBFS)
% function [pRev,Level,outputs] = demo_raa(filename,dBFS)
%
% 1. Description:
% 
%    Input parameters:
%       filename - EXPLAIN FILENAME
%       dBFS - EXPLAIN DBFS
%       insig - Input signal. Can have one channel (Left = Right, i.e., 
%               a diotic condition is used) or two channels.
% 
%    Output parameters:
%       pRev - reverberance estimates in Model Units
%       Level - array with overall levels dB:
%           Level(1) is Leq in dB(A)
%           Level(2) is Lmax in dB(A)
%           Level(3) is Leq in dB linear or dB(Z)
%           Level(4) is Lmax in dB linear or dB(Z)
%       outputs - struct containing additional information
%           outputs.pClar - provides the estimate of perceptual clarity
% 
% 2. Stand-alone example:
%       dBFS = 119; 
%       filename = '/home/alejandro/Documenten/Databases/dir01-Instruments/MBBM/lnr/Odeon_Orchestra/06-Bern_noModification/M127126-Odeon_mit-Musikern_ohne-Modifikation.ConvAural06.wav';
%       [pRev1,Level1] = demo_raa(filename,dBFS);
%
%       dBFS = 119;
%       filename = '/home/alejandro/Documenten/Databases/dir01-Instruments/MBBM/lnr/Odeon_Orchestra/08-Sydney_abs/160203_ohne-Musiker_nur-Donuts_02.ConvAural06.Wav';
%       [pRev2,Level2] = demo_raa(filename,dBFS);
% 
% 3. Additional info:
%       See reference Osses2017
%       Tested cross-platform: Yes
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2017
% Created on    : 21/01/2017
% Last update on: 21/01/2017 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    dBFS = dbspl(1); % 100 dB SPL attributed to an amplitude of 1, AMT default
end

Level = zeros(1,5); % Initialisation

[insig,fs] = audioread(filename);

subfs = 11025; % Hz

% insig = Create_sin(1000,1,fs);
% insig = setdbspl(insig,94,dBFS);
% insig = [insig insig];

insig = insig(1:round(10*fs),:);

insig = gaindb(1,dBFS-100)*insig; % gaindb replaces 'my' From_dB

dB = Do_SLM(insig,fs,'A','f',100);
dB = 10*log10(0.5*abs(10.^(dB(:,1)/10)+10.^(dB(:,2)/10)));
    
Level(1) = Get_Leq(dB); % LAeq
Level(2) = max(dB); % LAmax

dB = Do_SLM(insig,fs,'Z','f',100);
dB = 10*log10(0.5*abs(10.^(dB(:,1)/10)+10.^(dB(:,2)/10)));
Level(3) = Get_Leq(dB); % LZeq
Level(4) = max(dB); % LZmax

[outsig, fc, par] = dorp2011preproc(insig, fs, 'no_binaural','subfs',subfs);

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

disp('')

