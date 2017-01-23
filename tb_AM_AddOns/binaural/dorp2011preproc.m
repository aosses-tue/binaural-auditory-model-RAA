function [outsig, fc, par, Psi_dir, Psi_rev] = dorp2011preproc(insig, fs, varargin)
% function [outsig, fc, par, Psi_dir, Psi_rev] = dorp2011preproc(insig, fs, varargin)
%
% 1. Description:
%  Auditory model from van Dorp et. al. 2011
%   Usage: [outsig, fc, par] = dorp2011preproc(insig,fs);
%          [outsig, fc, par] = dorp2011preproc(insig,fs,...);
%
%   Input parameters:
%        insig  : input acoustic signal.
%        fs     : sampling rate.
%  
%   The input must have dimensions time x left/right channel x signal no.
%   The output has dimensions time x frequency x signal no. 
%  
%   Output parameters:
%       outsig  : output signal scaled in Model Units.
%       fc      : centre frequencies of the filter bank.
%       par     : parameters returned by the central processor:
%           - fc: centre frequencies of the filter bank.
%           - fs: sampling frequency of the auditory stream (downsampled to 
%                 16 kHz by default)
%           - numBands: number of audio frequency bands used by the model.
%           - numSamples: number of samples of the auditory stream.
%           - t_psi: time stamp of the auditory stream [s].
%           - PsiL: auditory stream of the  left channel (channel 1).
%           - PsiR: auditory stream of the right channel (channel 2).
%           - FL: global foreground level (direct stream).
%           - BL: global background level (reverberant stream).
%           - TL: total level = FL+BL.
%           - Llow: level of the low-frequency channels.
%           - TL_energy: relative energy per frequency band [%]
%           - numFrames: number of analysis frames.
%           - t_frame: time stamp of the frame-based estimates.
%           - FL_frame: foreground level per analysis frame.
%           - BL_frame: background level per analysis frame.
%           - TL_frame: total level per analysis frame.
%           - pRev
%           - pRev_frame
%           - sRev
%           - pClar
%           - pClar_frame
%           - sClar
%           - pASW
%           - sASW
%           - pLEV
%           - sLEV
%           - Psimin
%           - Psimin_dip
% 
%   The Dorp 2011 model consists of the following stages:
%   
%   1) Outer- and middle- ear filtering: band-pass filter with pass-band 
%      between 1000 and 4000 Hz. The slopes are approximately 6 dB/Oct.
%   2) a gammatone filter bank with 1-ERB spaced filters.
%   3,4) an envelope extraction stage done by half-wave rectification
%      followed by low-pass filtering to 770 Hz.
%   5) Absolute threshold of hearing: in this implementation this was omitted
%   6) an adaptation stage modelling nerve adaptation by a cascade of 5 loops.
%   7) Signal smoothing by applying a single-pole low-pass filter with a time
%      constant of 20-ms (fcut = 8 Hz; as in dau1996, breebaart2001). 
%   8) an excitation-inhibition (EI) cell model (not implemented yet).
%
%   Parameters for AUDITORYFILTERBANK, IHCENVELOPE, ADAPTLOOP and EICELL 
%   can be passed at the end of the line of input arguments.
%
%   Examples
%   --------
%
%   The following code sets up a simple test example :
%     % Setup parameters
%     fs      = 44100;            % Sampling rate
%     T       = 1;              % Duration
%     Spl1    = 75;               % SPL of input signal 1
%     Spl2    = 75;               % SPL of input signal 2
%     rho = 0;
% 
%     % Generate signals:
%     t  = [0:1/fs:T];
%     n1 = setdbspl(randn(length(t),1),Spl1);
%     n2 = setdbspl(randn(length(t),1),Spl2);
%     x1 = n1*sqrt((1+rho)/2) + n2*sqrt((1-rho)/2);
%     x2 = n1*sqrt((1+rho)/2) - n2*sqrt((1-rho)/2);
%
%     % % Run the model and plot it
%     [ei_map, fc] = dorp2011preproc([x1,x2], fs);
%     [idx]        = dorp2011preproc([x1,x2], fs,'central_envelopment');
% 
%     file = 'D:\Databases\dir01-Instruments\Piano\00-Original-files\pressionexpeCd5.wav';
%     [insig fs] = Wavread(file);
%     lvl = rmsdb(insig)+100;
%     [par psi] = raa(file,100);
% 
%   See also: headphonefilter_Dorp2011, auditoryfilterbank, ihcenvelope, adaptloop
%
%   References: breebaart2001binaural

%   AUTHOR : Alejandro Osses V., adapted partially from a code programmed
%            by Peter L. Soendergaard
  
% ------ Checking of input parameters ------------

if nargin<2
  error('%s: Too few input arguments.',upper(mfilename));
end;

if ~isnumeric(insig) 
  error('%s: insig must be numeric.',upper(mfilename));
end;

if ~isnumeric(fs) || ~isscalar(fs) || fs<=0
  error('%s: fs must be a positive scalar.',upper(mfilename));
end;

definput.import = {'dorp2011','auditoryfilterbank_','ihcenvelope','adaptloop_'};
definput.importdefaults={'gtf_dorp','ihc_breebaart','adt_dorp'};
definput.keyvals.subfs=[];

[flags,keyvals,flow,fhigh,basef]  = ltfatarghelper({'flow','fhigh','basef'},definput,varargin);

%% Truncates the signal to multiples of one second:
% This makes the global estimates equivalent to the mean of the framed-based
% values (if framelen is 5 s and hopsize 1 s, then the mean of the estimates
% 1, 6, 11... until the end (in steps of framelen/hopsize) provide this
% equivalence.
if flags.do_multiple_one_second
    Nactual = size(insig,1);
    N = floor(Nactual/fs)*fs; % N samples, rounded below to have only multiples of 1 s
    insig = insig(1:N,:);

    if Nactual-N ~= 0
        fprintf('The input signal has a length that is not a multiple of 1 s\n\t discarding %.0f samples\n',Nactual-N);
        if N == 0
            error('The input signal has to be at least 1 second long');
        end
    end
end

% ------ do the computation -------------------------
if flags.do_exclude
    
    dBFS = 100;
    dBA = Do_SLM(insig,fs,'A','f',dBFS);
    if size(insig,2) == 2
        dBA = 10*log10(0.5*abs(10.^(dBA(:,1)/10)+10.^(dBA(:,2)/10)));
        bDiotic = 0;
    else
        bDiotic = 1;
    end
    
    LAeq     = Get_Leq( dBA,fs);
    LAeq1sec = Get_Leq( dBA,fs, keyvals.hopsize); % 1-sec LAeq

    idx_no_silence = find( (LAeq1sec) >  LAeq-30 );
    
    keyvals.N              = length(LAeq1sec)      *keyvals.hopsize*fs;
    keyvals.N_no_silence   = length(idx_no_silence)*keyvals.hopsize*fs;
    keyvals.idx_no_silence = idx_no_silence;
else
    if size(insig,2) == 2
        bDiotic = 0;
    else
        bDiotic = 1;
    end
    
    keyvals.N              = N;
    keyvals.N_no_silence   = N;
    keyvals.idx_no_silence = 1:size(insig,1)/(keyvals.hopsize*fs);
end

if flags.do_outermiddleear
    % First-order Butterworth HP (cut-off = 1000 Hz) + LP (cut-off = 4000 Hz)
    [b1, b2, a1, a2] = headphonefilter_Dorp2011(fs,keyvals.order);
    outsig = filter(b1,a1,insig); % high-pass
    outsig = filter(b2,a2,outsig); % low-pass
    
elseif flags.do_nooutermiddleear
    outsig = insig;
    
end

%% Downsampling, if resquested:
% The downsampling at this stage, saves more computer resources (computations
% are done with less data):
if ~isempty(keyvals.subfs)
    outsig = fftresample(outsig,round(length(outsig)/fs*keyvals.subfs));
    fs = keyvals.subfs;

    if flags.do_exclude
        keyvals.N            = length(LAeq1sec)      *keyvals.hopsize*fs;
        keyvals.N_no_silence = length(idx_no_silence)*keyvals.hopsize*fs;
    else
        keyvals.N            = size(outsig,1);
        keyvals.N_no_silence = keyvals.N;
    end
end;

%% Apply the auditory filterbank
%       4th-order (default) complex filter bank. It matches the definition
%       presented by van Dorp. See Dorp2011, pag. 60 (out of 248): [...] Note 
%       that Breebaart proposed using third-order gammatone filters, but here a 
%       fourth-order filterbank is used [...]
%       In Breebaart2001, appendix 4.C explicit reference to the complex
%       gammatone filterbank is made.
[outsig, fc] = auditoryfilterbank(outsig,fs,'argimport',flags,keyvals);

%% 'haircell' envelope extraction
outsig = ihcenvelope(outsig,fs,'argimport',flags,keyvals);

%% non-linear adaptation loops
% keyvals.limit = 10;
outsig = adaptloop(outsig,fs,'argimport',flags,keyvals);

%% Binaural processor (Envelopment)
par = [];
if flags.do_binaural & bDiotic == 0
    % Not validated yet
end

%% 8-Hz Low-pass filter
[mlp_b, mlp_a, tc] = IRIfolp(8,fs); % 8 Hz LPF
outsig = filter(mlp_b,mlp_a,outsig);
   
%% Central processor:
% We set psi following RAA nomenclature
par.bDiotic    = bDiotic;
par.fc         = transpose(fc);
par.fs         = fs;
par.t_psi      = transpose(( 1:size(outsig,1) )/fs);

if bDiotic == 0
    par.PsiL       = outsig(:,:,1);
    par.PsiR       = outsig(:,:,2);
else
    par.PsiL       = outsig(:,:);
    par.PsiR       = outsig(:,:);
end
par.numBands   = size(par.PsiL,2);
par.numSamples = keyvals.N_no_silence; % it excludes silent segments if found

switch nargout
    case 3
        par = dorp2011centralprocessor(par,fs,flags,keyvals);
    case 4
        [par, Psi_dir] = dorp2011centralprocessor(par,fs,flags,keyvals);
    case 5 
        [par, Psi_dir, Psi_rev] = dorp2011centralprocessor(par,fs,flags,keyvals);
end

keyvals.idx_no_silence(keyvals.idx_no_silence > par.numFrames)=[];
par.IdxFrames = keyvals.idx_no_silence;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
