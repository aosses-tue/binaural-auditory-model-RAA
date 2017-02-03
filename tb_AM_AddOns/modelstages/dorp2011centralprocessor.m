function [par, Psi_dir, Psi_rev] = dorp2011centralprocessor(psi,fs,flags,keyvals)
% function [par, Psi_dir, Psi_rev] = dorp2011centralprocessor(psi,fs,flags,keyvals)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%       See also: Get_RAA_estimates.m (previous version of this file)
%       % TODO: frame-based estimates for ITD (basically do Llow and std frame-based)
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 12/07/2016
% Last update on: 21/07/2016 
% Last use on   : 21/07/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. Fixed parameters:
Psimin     =  keyvals.Psimin;  
Tmin       =  keyvals.Tmin;
Psimin_dip = -keyvals.Psimin_dip; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bDiotic = psi.bDiotic;
PsiL = psi.PsiL;
if bDiotic == 0
    % Normal stereo configuration:
    PsiR = psi.PsiR; 
else
    % if bDiotic == 1, then Left is copied to Right channel:
    PsiR = PsiL;
end
    
N = keyvals.N_no_silence; % coming from keyvals.N

[PsiL_dir, PsiL_rev,LpsiL] = il_split_stream(PsiL,fs,N,Psimin,Psimin_dip,Tmin);
if bDiotic == 0
    [PsiR_dir, PsiR_rev,LpsiR] = il_split_stream(PsiR,fs,N,Psimin,Psimin_dip,Tmin);
else
    PsiR_dir = PsiL_dir;
    PsiR_rev = PsiL_rev;
    LpsiR = LpsiL;
end

[TL_frame,t_fr,numFr,TL,TL_ene] = il_get_frame_calc(PsiL,PsiR,fs,keyvals); % Total level
[FL_frame,t_fr,numFr,FL,FL_ene] = il_get_frame_calc(PsiL_dir,PsiR_dir,fs,keyvals); % Foreground level
[BL_frame,t_fr,numFr,BL,BL_ene] = il_get_frame_calc(PsiL_rev,PsiR_rev,fs,keyvals); % Background level
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargout > 1
    % Clears memory if possible
    Psi_dir(:,:,1) = PsiL_dir; clear PsiL_dir;
    Psi_dir(:,:,2) = PsiR_dir; clear PsiR_dir;
end
if nargout > 2
    % Clears memory if possible
    Psi_rev(:,:,1) = PsiL_rev; clear PsiL_rev;
    Psi_rev(:,:,2) = PsiR_rev; clear PsiR_rev;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 4.1 Reverberance
pREV = BL;

%%% 4.2 Clarity:
pCLA = FL/BL;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

par.fc         = psi.fc;
par.fs         = psi.fs; 
par.t_psi      = psi.t_psi;
par.numSamples = psi.numSamples;
par.numBands   = psi.numBands;
par.numFrames  = numFr;
par.FL         = FL;
par.BL         = BL;
par.TL         = TL;
par.TL_energy  = TL_ene; % percentage of energy per auditory filter
par.BL_energy  = BL_ene; 
par.FL_energy  = FL_ene; 
par.FL_frame   = FL_frame;
par.BL_frame   = BL_frame;
par.TL_frame   = TL_frame;
par.t_frame    = t_fr;

par.pRev       = pREV;
par.pRev_frame = BL_frame; % as suggested by van Dorp
par.pClar      = pCLA;
par.pClar_frame= FL_frame./BL_frame; % as suggested by van Dorp
par.Psimin     = Psimin;
par.Psimin_dip = Psimin_dip;
par.LpsiL      = LpsiL;
par.LpsiR      = LpsiR;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EOF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inline functions:
function [outDir,outRev,Lpsi] = il_split_stream(Stream,fs,N,mumin,mumin_dip,Tmin)

outDir  = zeros(size(Stream));
outRev  = zeros(size(Stream));

Lpsi = 1/N*sum(abs(Stream)); % Eq. 3.15

for i = 1:size(Stream,2)
    
    Psimin     =     mumin*Lpsi(i); % Eq. 3.14
    Psimin_dip = mumin_dip*Lpsi(i);
    
    idx = find(Stream(:,i) >= Psimin);
    idxDir_above = il_detect_segment(idx,fs,Tmin);
    idx = find(Stream(:,i) < Psimin_dip);
    
    idxDir_below = il_detect_segment(idx,fs,Tmin);
    
    idxDir = sort([idxDir_above idxDir_below]);
        
    idxRev = 1:size(Stream(:,i),1);
    idxRev(idxDir) = [];
    
    outDir(idxDir,i) = Stream(idxDir,i);
    outRev(idxRev,i) = Stream(idxRev,i);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function idxDir = il_detect_segment(idx,fs,Tmin)

[iblocks(idx), idxL, idxU] = il_detect_blocks(idx);

try
    idxL = idx(idxL); % convert idxL relative to 'iblocks'
    idxU = idx(idxU);
    dur = (idxU-idxL)/fs;
    iBlocks2use = find(dur(:) >= Tmin);
catch
    iBlocks2use = [];
    % disp('No blocks found at least in one frequency band...')
end

idxDir = [];
for k = 1:length(iBlocks2use)
    idx2usetmp = find(iblocks==iBlocks2use(k));
    idxDir = [idxDir idx2usetmp];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [iblocks, idxL, idxU] = il_detect_blocks(idx)

nConsecutive  =diff(idx);
iblocks = [];

idxBound = find(nConsecutive>1);
idxBound = [idxBound; length(idx)];

if idxBound(1) ~= 1
    idxBound = [1; idxBound];
end

for i = 1:length(idxBound)-1
    idxL(i) = idxBound(i);
    idxU(i) = idxBound(i+1)-1;
    iblocks(idxL(i):idxU(i)) = i;
end

if length(iblocks) ~=0
    iblocks(end+1) = iblocks(end);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [out, t, NumFrames, Level, Level_energy] = il_get_frame_calc(PsiL,PsiR,fs,keyvals)

Stream  = sqrt(PsiL.^2+ PsiR.^2);
K       = size(Stream,2);
N_no_silence = keyvals.N_no_silence;
N       = keyvals.N;

Level        =       1/(N_no_silence*K)*sum(sum(Stream));
Level_energy = 100*( 1/(N_no_silence*K)*sum(Stream) )/Level;

framelenS = round(keyvals.framelen*fs);
hopsizeS  = round(keyvals.hopsize *fs);
NumFrames = floor((N-framelenS)/hopsizeS)+1;

if N_no_silence-framelenS < 0
    error('Signal is too short: signal is shorter than the analysis frame length');
end
    
for i = 1:NumFrames
    idxi = (i-1)*hopsizeS+1;
    idxf = idxi+framelenS-1;
    t(i) = (idxi-1)*fs; % t(1) is assigned to 0 s
    out(i) = 1/(framelenS*K)*sum(sum(Stream(idxi:idxf,:)));
end
