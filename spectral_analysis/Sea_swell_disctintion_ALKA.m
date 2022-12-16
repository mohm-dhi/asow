% Sea and swell disctintion
% ALKA - Nov 2022
clear all; close all; clc

%% SETTINGS
MainDir = '\\usden1-nas3\41806379\';
Outdir = 'D:\PTFI';
loc = 'Pipeline';

%% LOAD M_STRUCTURES
disp('Loading Spectra');
load([MainDir,'06_Results\M_STRUCTURES\',loc,'_2DSpectra.mat'])

disp('Loading Wind');
load([MainDir,'06_Results\M_STRUCTURES\Offshore_WD.mat']) % wind
load([MainDir,'06_Results\M_STRUCTURES\Offshore_WS.mat'])

%% ENSURE OVERLAP IN TIME
[time,Ia,Ib] = intersect(WS.time, ED2f.time);

WS.time = WS.time(Ia);
WS.data = WS.data(Ia);
WS.datetime = WS.datetime(Ia);
WS.ttt = [time(1), time(end), 60 60];

WD.time = WD.time(Ia);
WD.data = WD.data(Ia);
WD.datetime = WD.datetime(Ia);
WD.ttt = [time(1), time(end), 60 60];

ED2f.time = ED2f.time(Ib);
ED2f.data = ED2f.data(Ib,:,:);
ED2f.datetime = ED2f.datetime(Ib);
ED2f.ttt = [time(1), time(end), 60 60];

% define waterlevels
WL = 0; %

%% RUN PARTITION
cd(Outdir)
m_waveage(ED2f,WL,WS,WD,'factor',0.78,'power',0.2);

%% CALCULATE PARAMETERS
load([loc,'_Waveage_ED2f_deg_SW.mat'])

labels = {'total','sea','swell'};

for k = 1:length(labels)
    
    ED2f_Part = M;
    ED2f_Part.name = ['Pipeline_',labels{k}];
    ED2f_Part.item = ED2f_Part.item{k};
    ED2f_Part.label = ED2f_Part.label{k};
    ED2f_Part.data = ED2f_Part.data(:,:,:,k);
    
    % m_spectra(ED2f_Part,'plotflag',0)
    eval([ED2f_Part.name,'= m_2DspecInt(ED2f_Part,ED2f_Part.xaxis,ED2f_Part.yaxis)']);
    save([ED2f_Part.name,'.mat'],ED2f_Part.name)
end

return
%% RUN WATERSHED
[AA,Ef]=partition(freq,dir,E,wfc,fw,sw);

%%
% fig = figure;
% subplot(1,3,1)
% plot(FSRU_total.Values(:,1),FSRU_total.Values(:,2),'.')
% xlim([0,5]);ylim([0,15])
% title('total'); xlabel('Hs'); ylabel('Tp')
% subplot(1,3,2)
% plot(FSRU_sea.Values(:,1),FSRU_sea.Values(:,2),'.')
% xlim([0,5]);ylim([0,15])
% title('sea'); xlabel('Hs'); ylabel('Tp')
% subplot(1,3,3)
% plot(FSRU_swell.Values(:,1),FSRU_swell.Values(:,2),'.')
% xlim([0,5]);ylim([0,15])
% title('swell'); xlabel('Hs'); ylabel('Tp')

% subplot(2,3,4)
% plot(FSRU_total.Values(:,7),FSRU_total.Values(:,2),'.')
% xlim([150,360]);ylim([0,15])
% title('total'); xlabel('MWD'); ylabel('Tp')
% subplot(2,3,5)
% plot(FSRU_sea.Values(:,7),FSRU_sea.Values(:,2),'.')
% xlim([150,360]);ylim([0,15])
% title('sea'); xlabel('MWD'); ylabel('Tp')
% subplot(2,3,6)
% plot(FSRU_swell.Values(:,7),FSRU_swell.Values(:,2),'.')
% xlim([150,360]);ylim([0,15])
% title('swell'); xlabel('MWD'); ylabel('Tp')