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
  
f = mfilename('fullpath');
flast = il_strsplit(f,il_delim);
L = length(flast{end});

bp = f(1:end-L);
disp('')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = il_strsplit(s,d)

str={};
prev=1;
count=1;
for i=1:length(s)
    if (s(i)==d)
        if (prev<=i-1)
            str{count}=s(prev:i-1);
        end
        count=count+1;
        prev=i+1;
    end
end

str{count}=s(prev:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function delimiter = il_delim

if isunix
    delimiter = '/';
else
    delimiter = '\';
end
