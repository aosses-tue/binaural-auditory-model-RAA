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
if bDiotic
    PsiR = psi.PsiR; 
else
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
% Low-frequency energy (not used for reverberation and clarity estimates).

if flags.do_binaural
    itdFrameLen = 50e-3; % ms
    [ITD_dir,ITD_rev] = il_split_ITD_stream(psi.t_itd,psi.itd,psi.t_psi,PsiL_dir,PsiR_dir,itdFrameLen);
    sigma_dir = std(ITD_dir);
    sigma_rev = std(ITD_rev);

    muASW = keyvals.muASW;
    nuASW = keyvals.nuASW;
    muLEV = keyvals.muLEV;
    nuLEV = keyvals.nuLEV;
    q0 = 9;
    q1 = 20;
    idx_i = q0-4;
    idx_f = q1-4;
    Q = (idx_f-idx_i)+1;
    sigma_dir = 1/Q*sum(sigma_dir(idx_i:idx_f)); % Dorp2011, Eq. 3.24
    sigma_rev = 1/Q*sum(sigma_rev(idx_i:idx_f));

    z0    = 5;
    z1    = 9;
    idx_i = z0-4; % band number (initial)
    idx_f = z1-4; % band number (final)
    idx4Llow = idx_i:idx_f;
    Kl    = (idx_f-idx_i)+1; % resulting number of frequency bands, for the low-freq calculation

    % Dorp2011, Eq. 3.25, pp 82:
    Llow = 1/(Kl*N)*sum(  sum( sqrt((PsiL(:,idx4Llow)).^2+(PsiR(:,idx4Llow)).^2 ))  );
    clear PsiL PsiR; % to free memory

    pASW = muASW*Llow+log10( 1+nuASW*sigma_dir*1e3 ); % Dorp2011, Eq. 3.23
    pLEV = muLEV*BL  +log10( 1+nuLEV*sigma_rev*1e3 ); % Dorp2011, Eq. 3.26
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargout > 1
    Psi_dir(:,:,1) = PsiL_dir; clear PsiL_dir;
    Psi_dir(:,:,2) = PsiR_dir; clear PsiR_dir;
end
if nargout > 2
    Psi_rev(:,:,1) = PsiL_rev; clear PsiL_rev;
    Psi_rev(:,:,2) = PsiR_rev; clear PsiR_rev;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    rev   cla   asw  lev
a = [0.207 1.17  2.76 2.41];
b = [21.4  2.20  2.84 2.64];

%%% 4.1 Reverberance
pREV = BL;
sREV = 1/(1+exp(-a(1)*(pREV-b(1))));

%%% 4.2 Clarity:
pCLA = FL/BL;
sCLA = 1/(1+exp(-a(2)*(pCLA-b(2))));

if flags.do_binaural
    %%% 4.3 ASW
    sASW = 1/(1+exp(-a(3)*(pASW-b(3))));

    %%% 4.4 Listener Envelopment
    sLEV = 1/(1+exp(-a(4)*(pLEV-b(4))));
end

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
par.sRev       = sREV;
par.pClar      = pCLA;
par.sClar      = sCLA;
if flags.do_binaural
    par.Llow   = Llow;
    par.pASW   = pASW;
    par.sASW   = sASW;
    par.pLEV   = pLEV;
    par.sLEV   = sLEV;
end
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
    try
        idxDir_below = il_detect_segment(idx,fs,Tmin);
    catch
        disp('')
    end
    idxDir = sort([idxDir_above idxDir_below]);
        
    idxRev = 1:size(Stream(:,i),1);
    idxRev(idxDir) = [];
    
    outDir(idxDir,i) = Stream(idxDir,i);
    outRev(idxRev,i) = Stream(idxRev,i);
    
    disp('')
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
% tic
for i = 1:length(idxBound)-1
    idxL(i) = idxBound(i);
    idxU(i) = idxBound(i+1)-1;
    iblocks(idxL(i):idxU(i)) = i;
end
% toc

try
    if length(iblocks) ~=0
        iblocks(end+1) = iblocks(end);
    end
catch
    warning('Check block detection code')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ITD_dir,ITD_rev] = il_split_ITD_stream(t_itd,itd,t_psi,PsiL_dir,PsiR_dir,itdFrameLen)

idxL_dir = (PsiL_dir ~= 0);
idxR_dir = (PsiR_dir ~= 0);
zeta_DIR_t = idxL_dir+idxR_dir;

t_idx(1,1) = 1;

for i = 1:size(t_itd,1)
    if i ~= 1
        t_idx(i,1) = find(t_psi > t_itd(i) - itdFrameLen/2,1,'first')-1;
    end
    if i ~= size(t_itd,1)
        t_idx(i,2) = find(t_psi < t_itd(i) + itdFrameLen/2,1,'last')-1;
    else
        t_idx(i,2) = size(t_psi,1);
    end
    
    % j = 5;
    zeta_DIR(i,:) = 0.5*round(mean(zeta_DIR_t(t_idx(i,1):t_idx(i,2),:)));
    disp('')
    
end
ITD_dir =    zeta_DIR .*itd; % Eq 3.17
ITD_rev = (1-zeta_DIR).*itd; % Eq 3.18

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [out, t, NumFrames, Level, Level_energy] = il_get_frame_calc(PsiL,PsiR,fs,keyvals)

Stream  = sqrt(PsiL.^2+ PsiR.^2);
K       = size(Stream,2);
N_no_silence = keyvals.N_no_silence;
N       = keyvals.N;

Level        = 1/(N_no_silence*K)*sum(sum(Stream));
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
