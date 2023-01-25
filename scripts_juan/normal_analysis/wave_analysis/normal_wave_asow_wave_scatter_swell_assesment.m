%%
clear all; clc

addpath(genpath('C:\DHI\Projects\AtlanticShores\Scripts\asow\scripts_juan\normal_analysis'));

%%

fdir = 'C:\DHI\Projects\AtlanticShores\Data\TimeSeries\Deliverable\CSVFinal\structs\';
odir_scatter = 'C:\DHI\Projects\AtlanticShores\Analysis\Normal\Waves\Results\SwellAssesment\';
main_dir = 'Swell_Assesment';
no_locs = 7;

%% Hm0Sea/Hm0Total Omni

params1 = {'Hm0_Total'};
params2 = {'Hm0_Sea'};
% direc = {'PWD_Total','PWD_Sea','PWD_Swell'};
bins_hm0_total = {0:2:8, 0:2:8, 0:2:8};
bins_hm0_sea = {0:2:10, 0:5:15, 0:5:15};
% leg = {'H_{m0,Total}','H_{m0,Sea}','H_{m0,Swell}'};

for i=1:no_locs

    cd([odir_scatter 'P' num2str(i)]);
    mkdir('Sea_Total');
    cd('Sea_Total');

    % load structure
    fname = [fdir 'ASOW' num2str(i) '_all_structs.mat'];
    load(fname);

    hm0_total = asow_params.Hm0_Total;
    hm0_total.data(hm0_total.data<0)=NaN;
    hm0_sea = asow_params.Hm0_Sea;
    hm0_sea.data(hm0_sea.data<0)=NaN;

    % modify struct items for plotting
    hm0_total.bins = bins_hm0_total{1};
    hm0_sea.bins = bins_hm0_sea{1};
    hm0_total.legend = 'SW_{US-EC}';
    hm0_sea.legend = 'SW_{US-EC}';
    hm0_total.ttt_str_long = ' (1979-01-15–2021-12-31; \Deltat=1h) ';
    hm0_total.ttt = [datenum('1979-01-15') datenum('2021-12-31') 60];
    hm0_total.xyz_str =  ['(' num2str(-1*hm0_total.xyz(1),'%.3f') 'W; ' num2str(hm0_total.xyz(2),'%.3f') 'N; ' num2str(hm0_total.xyz(3),'%.1f') 'mMSL' ')'];

    % do statistics
    %         [fitvals] = m_scatter(hm0_struct,tm02_struct,'density','quantiles',[0.05 0.5 0.95],'fit_func_Q','Poly','nolegend');
    %         save('fitvals.mat','fitvals');
    [fitvals] = m_scatter(hm0_total,hm0_sea,'plottype','scatter','density','quantiles',[0.05 0.5 0.95],'fit_func_Q','Poly');
    save('fitvals.mat','fitvals');
    cd('../..');

end

%% Hm0Swell/Hm0Sea Omni

bins_hm0_swell = {0:2:8, 0:2:8, 0:2:8};
bins_hm0_sea = {0:2:10, 0:5:15, 0:5:15};

for i=1:no_locs

    cd([odir_scatter 'P' num2str(i)]);
    mkdir('Swell_Sea');
    cd('Swell_Sea');

    % load structure
    fname = [fdir 'ASOW' num2str(i) '_all_structs.mat'];
    load(fname);

    hm0_sea = asow_params.Hm0_Sea;
    hm0_sea.data(hm0_sea.data<0)=NaN;
    hm0_swell = asow_params.Hm0_Swell;
    hm0_swell.data(hm0_swell.data<0)=NaN;

    % modify struct items for plotting
    hm0_sea.bins = bins_hm0_swell{1};
    hm0_swell.bins = bins_hm0_sea{1};
    hm0_sea.legend = 'SW_{US-EC}';
    hm0_swell.legend = 'SW_{US-EC}';
    hm0_sea.ttt_str_long = ' (1979-01-15–2021-12-31; \Deltat=1h) ';
    hm0_sea.ttt = [datenum('1979-01-15') datenum('2021-12-31') 60];
    hm0_sea.xyz_str =  ['(' num2str(-1*hm0_sea.xyz(1),'%.3f') 'W; ' num2str(hm0_sea.xyz(2),'%.3f') 'N; ' num2str(hm0_sea.xyz(3),'%.1f') 'mMSL' ')'];

    % do statistics
    %         [fitvals] = m_scatter(hm0_struct,tm02_struct,'density','quantiles',[0.05 0.5 0.95],'fit_func_Q','Poly','nolegend');
    %         save('fitvals.mat','fitvals');
    [fitvals] = m_scatter(hm0_sea,hm0_swell,'plottype','scatter','density','quantiles',[0.05 0.5 0.95],'fit_func_Q','Poly');
    save('fitvals.mat','fitvals');
    cd('../..');

end

%% Hm0Swell/Hm0Total Omni

bins_hm0_swell = {0:2:8, 0:2:8, 0:2:8};
bins_hm0_sea = {0:2:10, 0:5:15, 0:5:15};

for i=1:no_locs

    cd([odir_scatter 'P' num2str(i)]);
    mkdir('Swell_Total');
    cd('Swell_Total');

    % load structure
    fname = [fdir 'ASOW' num2str(i) '_all_structs.mat'];
    load(fname);

    hm0_total = asow_params.Hm0_Total;
    hm0_total.data(hm0_total.data<0)=NaN;
    hm0_swell = asow_params.Hm0_Swell;
    hm0_swell.data(hm0_swell.data<0)=NaN;

    % modify struct items for plotting
    hm0_total.bins = bins_hm0_swell{1};
    hm0_swell.bins = bins_hm0_sea{1};
    hm0_total.legend = 'SW_{US-EC}';
    hm0_swell.legend = 'SW_{US-EC}';
    hm0_total.ttt_str_long = ' (1979-01-15–2021-12-31; \Deltat=1h) ';
    hm0_total.ttt = [datenum('1979-01-15') datenum('2021-12-31') 60];
    hm0_total.xyz_str =  ['(' num2str(-1*hm0_total.xyz(1),'%.3f') 'W; ' num2str(hm0_total.xyz(2),'%.3f') 'N; ' num2str(hm0_total.xyz(3),'%.1f') 'mMSL' ')'];

    % do statistics
    %         [fitvals] = m_scatter(hm0_struct,tm02_struct,'density','quantiles',[0.05 0.5 0.95],'fit_func_Q','Poly','nolegend');
    %         save('fitvals.mat','fitvals');
    [fitvals] = m_scatter(hm0_total,hm0_swell,'plottype','scatter','density','quantiles',[0.05 0.5 0.95],'fit_func_Q','Poly');
    save('fitvals.mat','fitvals');
    cd('../..');

end
