function exp_osses2017
% function exp_osses2017
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
%           demo_osses2017; % ... and follow the instructions on the screen
% 
% 3. Additional info:
%       See reference Osses2017
%       Tested cross-platform: Yes
% 
% Old name in RAA v1.0: demo_osses2017.m
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2017
% Created on    : 21/01/2017
% Last update on: 03/02/2017 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compatible version to be uploaded to GIT repository
bCalc = input([mfilename ':\n - Type 0 if you want to plot the stored results (choose this option if you want to generate a resulting figure only), or\n - Type 1 if you want to run all the RAA simulations : ']);
bRead = ~bCalc; 

if bCalc
    
    dir_where = raabasepath; % Example path in Windows:  dir_where = 'D:\MATLAB_RAA\tb_AM_AddOns\';
                             % Example path in Unix sys: dir_where = '/home/alejandro/Documenten/MATLAB_RAA/tb_AM_AddOns/';
    dir_where = [dir_where 'tb_AM_AddOns' filesep 'auxdata' filesep 'osses2017' filesep];

    files = il_get_filenames;

    dBFS = 100; % 100 dB SPL attributed to an amplitude of 1, old AMT default

    N = length(files);
    for i = 1:N
        fprintf('Processing file %.0f of %.0f\n',i,N);
        [pRev(i,:),Level(i,:),outputs] = osses2017([dir_where files{i}],dBFS);
    end

    pRev_est = reshape(pRev(:,1),8,6);
    
    disp('The following values of pRev_est should be exactly the same ones')
    disp(['as in ' mfilename ', lines 61-68 of the code: '])
    pRev_est
    
end

if bRead
    pRev_est = [     9.3727   14.7613   17.4746    6.2734    7.0316   17.3556; ...
                     8.1550   14.8859   14.8696    4.6825    6.6977   18.2143; ...
                     8.4383   14.6738   17.6120    5.6642    8.1008   19.3884; ...
                     8.7549   16.3137   17.6127    3.9633    6.6082   18.2909; ...
                     9.3461   16.8591   17.9122    5.7165    8.0234   18.6240; ...
                     9.8489   17.3096   18.7817    4.9929    7.4803   20.1715; ...
                    10.3667   18.1947   18.3527    5.4554    7.6701   21.5292; ...
                    13.8343   21.3438   23.3105    8.0382    9.6815   25.1301];
end

instr = {'Vio','Cello','DBass','Flute','Trum','Ti'};
[T30toti,EDTtoti] = il_get_Reverb;
pRev_paper        = il_get_pRev_paper;

try
    for i = 1:6
        [rpT(i,1),rpT(i,2)] = corr(pRev_est(:,i),T30toti(i,:)','type','Pearson');
        [rpE(i,1),rpE(i,2)] = corr(pRev_est(:,i),EDTtoti(i,:)','type','Pearson');
    end
catch
    warning('You don''t seem to have the statistical toolbox...')
end

legend4this = [];
offsetx = 0.1;
figure;
for i = 1:6
    x = [i-4*offsetx i-3*offsetx i-2*offsetx i-1*offsetx i i+1*offsetx i+2*offsetx i+3*offsetx];
    h1 = plot(x,pRev_paper(:,i),'ks-','MarkerFaceColor','k'); hold on
    h2 = plot(x,pRev_est(:,i),'ro-','MarkerFaceColor','r','LineWidth',2);
    legend4this{end+1} = instr{i};
end
set(gca,'XTick',1:6);
set(gca,'XTickLabel',legend4this);
grid on;
legend([h2 h1],'pRev as published','pRev from 10-s excerpts','Location','NorthWest');
ylim([2 28])
set(gca,'YTick',4:2:26)
ylabel('pRev [MU]')
h = gcf;
Pos = get(h,'Position');
Pos(3) = 1.6*Pos(3);
set(h,'Position',Pos);

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Musical instruments (rows 1 to 6 in the next data): ')
transpose(instr)
disp('(Each musical instrument had a different binaural room impulse response)')
disp('%%%')
disp('   ')
disp('EDT for Rooms 1 to 8 (columns 1 to 8): ')
EDTtoti

disp('T30 for Rooms 1 to 8 (columns 1 to 8): ')
T30toti
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(['EOF: ' mfilename])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function files = il_get_filenames
% function files = il_get_filenames
%
% File names of the 10-s excerpts uploaded to GITHUB

files = {'b01-Vio_1-room-1-Aabs.wav'; ...
         'b01-Vio_1-room-2-Babs.wav'; ...
        'b01-Vio_1-room-3-A.wav'; ...
        'b01-Vio_1-room-4-Cabs1.wav'; ...
        'b01-Vio_1-room-5-B.wav'; ...
        'b01-Vio_1-room-6-Cabs2.wav'; ...
        'b01-Vio_1-room-7-D.wav'; ...
        'b01-Vio_1-room-8-C.wav'; ...
        'b06-Cello-room-1-Aabs.wav'; ...
        'b06-Cello-room-2-Babs.wav'; ...
        'b06-Cello-room-3-A.wav'; ...
        'b06-Cello-room-4-Cabs1.wav'; ...
        'b06-Cello-room-5-B.wav'; ...
        'b06-Cello-room-6-Cabs2.wav'; ...
        'b06-Cello-room-7-D.wav'; ...
        'b06-Cello-room-8-C.wav'; ...
        'b07-DBass-room-1-Aabs.wav'; ...
        'b07-DBass-room-2-Babs.wav'; ...
        'b07-DBass-room-3-A.wav'; ...
        'b07-DBass-room-4-Cabs1.wav'; ...
        'b07-DBass-room-5-B.wav'; ...
        'b07-DBass-room-6-Cabs2.wav'; ...
        'b07-DBass-room-7-D.wav'; ...
        'b07-DBass-room-8-C.wav'; ...
        'b08-Flute-room-1-Aabs.wav'; ...
        'b08-Flute-room-2-Babs.wav'; ...
        'b08-Flute-room-3-A.wav'; ...
        'b08-Flute-room-4-Cabs1.wav'; ...
        'b08-Flute-room-5-B.wav'; ...
        'b08-Flute-room-6-Cabs2.wav'; ...
        'b08-Flute-room-7-D.wav'; ...
        'b08-Flute-room-8-C.wav'; ...
        'b20-Trum_1-room-1-Aabs.wav'; ...
        'b20-Trum_1-room-2-Babs.wav'; ...
        'b20-Trum_1-room-3-A.wav'; ...
        'b20-Trum_1-room-4-Cabs1.wav'; ...
        'b20-Trum_1-room-5-B.wav'; ...
        'b20-Trum_1-room-6-Cabs2.wav'; ...
        'b20-Trum_1-room-7-D.wav'; ...
        'b20-Trum_1-room-8-C.wav'; ...
        'b22-Ti-room-1-Aabs.wav'; ...
        'b22-Ti-room-2-Babs.wav'; ...
        'b22-Ti-room-3-A.wav'; ...
        'b22-Ti-room-4-Cabs1.wav'; ...
        'b22-Ti-room-5-B.wav'; ...
        'b22-Ti-room-6-Cabs2.wav'; ...
        'b22-Ti-room-7-D.wav'; ...
        'b22-Ti-room-8-C.wav'};

function [T30toti,EDTtoti] = il_get_Reverb
% function [T30toti,EDTtoti] = il_get_Reverb
%
% 1. Description:
%       Inline function with the reverberation values (T30 and EDT) derived
%       from the binaural impulse responses (BRIR) used to auralise each of the 
%       23 instrument groups of the Odeon Orchestra. The BRIRs were all different
%       since they account also for the directionality of each instrument
%       at the listener position.
%
%       T30toti and EDTtoti have 8 columns, where each column corresponds to
%       one acoustic condition in the following order: 
%           1. Room 1 - Aabs
%           2. Room 2 - Babs
%           3. Room 3 - A
%           4. Room 4 - Cabs1
%           5. Room 5 - B
%           6. Room 6 - Cabs2
%           7. Room 7 - D
%           8. Room 8 - C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T30toti = [ 1.3117    1.1310    1.8920    1.3135    1.9645    1.4445    2.1972    2.5072; ... % vio1a
            1.0017    1.1080    1.3315    1.0708    1.9457    1.3895    2.2740    2.5917; ... % cello
            1.1495    1.2360    1.2080    1.1578    2.0543    1.2897    2.1777    2.4855; ... % dbass
            1.8835    1.1777    1.2590    1.0430    2.0000    1.3015    2.2948    2.5855; ... % flute
            0.9307    1.2542    1.6060    1.0300    2.0250    1.2910    2.1310    2.5610; ... % trum1
            1.0797    1.2750    1.7083    1.0328    2.0118    1.2317    2.1270    2.4785];    % ti

EDTtoti =[  0.6368    1.0410    0.6627    0.9940    1.4900    1.2650    1.7062    2.4013; ... % vio1a
            0.7170    1.1285    0.7318    1.1115    1.3820    1.2175    1.5765    2.4653; ... % cello
            0.8680    0.8355    0.6658    0.9968    1.3670    1.2902    1.7695    2.4840; ... % dbass
            1.4635    0.8435    0.5367    1.0355    1.2268    1.2098    1.5722    2.4620; ... % flute
            0.5282    0.8835    0.5325    1.1250    1.4703    1.3530    1.3348    2.3965; ... % trum1
            0.4800    0.4765    1.2492    1.2020    0.8798    1.4750    1.8638    2.5030];    % ti

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Rev_ans = il_get_pRev_paper

Rev_ans = [ 12.2633   11.2291   12.0586   11.9042   13.0040   13.1479   12.9816   16.4841; ... % vio (vio1a, vio1b, vio2a, vio2b)
            13.6059   13.2556   14.3517   14.9818   15.8884   15.9189   16.3784   18.5668; ... % cello
            15.0666   12.2480   14.4881   15.0623   14.6228   16.1827   14.3101   20.7709; ... % dbass
            11.1275    8.5834    9.6291    7.8576    8.5442    8.9512    7.2913   11.8134; ... % flute
             6.4068    5.2814    5.8944    5.4206    5.0831    6.4100    5.0751    7.4184; ... % trum (trum1, trum2)
            18.2858   19.6861   20.8324   19.7431   20.8737   20.8191   21.5170   24.8661];    % ti

Rev_ans = transpose(Rev_ans);