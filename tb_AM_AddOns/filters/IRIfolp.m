function [b,a,timeconstant] = IRIfolp(f0,fs)
% function [b,a,timeconstant] = IRIfolp(f0,fs)
%
% 1. Description:
%       IRIfolp.m - digital first-order lowpass-filter transfer function
%                   (impulse response invariance)
%
%       Usage: [b,a] = folp(f0,fs)
%
%       f0 = cutoff frequency of the lowpass filter in Hz
%       fs = sampling rate in Hz
%
%       b = [b0,b1, ... ,bN] = numerator coefficients
%       a = [ 1,a1, ... ,aN] = denominator coefficients
%
% 2. Stand-alone example:
%       % If a time constant is required:
%       timeconstant = 20e-3; % 20 ms
%       f0 = 1/(2*pi*timeconstant);
%       fs = 44100;
%       [b a] = IRIfolp(f0,fs);
%       figure;
%       freqz(b,a,4096);
% 
%       % Almost the same example, but just giving f0:
%       f0 = 8; % Approx. frequency related to a time constant of 20 ms
%       fs = 44100;
%       [b a] = IRIfolp(f0,fs);
%       figure;
%       freqz(b,a,4096);
%        
% 3. Additional info:
%       Tested cross-platform: Yes
% 
% Obtained from Two!Ears team in 2015
% Comments by Alejandro Osses, HTI, TU/e, the Netherlands 2014-2015
% Last edit on: 16/08/2015
% Last use on : 16/08/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = exp( -(2*pi*f0)/fs );
b = 1 - a;
a = [1, -a];

timeconstant = 1/(2*pi*f0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eof
