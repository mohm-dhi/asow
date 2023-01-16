%%
clear all; clc

addpath(genpath('C:\DHI\Projects\AtlanticShores\Scripts\asow\scripts_juan\normal_analysis'));

%%

fdir = 'C:\DHI\Projects\AtlanticShores\Data\TimeSeries\Deliverable\CSVFinal\structs\';
odir_hmo_wshub = 'C:\DHI\Projects\AtlanticShores\Analysis\Normal\Waves\Scatter_Hm0_WSHub\';
no_locs = 7;

%% Hm0/Tp Omni

params1 = {'Hm0_Total','Hm0_Sea','Hm0_Swell'};
params2 = {'WindSpd_Hub'};
% direc = {'PWD_Total','PWD_Sea','PWD_Swell'};
bins_hm0 = {0:2:8, 0:2:8, 0:2:8};
bins_wshub = {0:5:45};

cd([odir_hmo_wshub 'Omni']);

for i=1:no_locs

    % load structure
    fname = [fdir 'ASOW' num2str(i) '_all_structs.mat'];
    load(fname);

    for p = 1:length(params1)

        mkdir(params1{p});
        cd(params1{p});
        mkdir(['ASOW' num2str(i)]);
        cd(['ASOW' num2str(i)])
        
        hm0_struct = asow_params.(params1{p});
        hm0_struct.data(hm0_struct.data<0)=NaN;
        wshub_struct = asow_params.(params2{1});
        wshub_struct.data(wshub_struct.data<0)=NaN;

        % modify struct items for plotting
        hm0_struct.bins = bins_hm0{p};
        wshub_struct.bins = bins_wshub{1};
        hm0_struct.legend = 'SW_{US-EC}';
        wshub_struct.legend = 'SW_{US-EC}';
        wshub_struct.label = 'WS_{Hub}';
        wshub_struct.unit = 'm/s';
        hm0_struct.ttt_str_long = ' (1979-01-15–2021-12-31; \Deltat=1h) ';
        hm0_struct.ttt = [datenum('1979-01-15') datenum('2021-12-31') 60];
        hm0_struct.xyz_str =  ['(' num2str(-1*hm0_struct.xyz(1),'%.3f') 'W; ' num2str(hm0_struct.xyz(2),'%.3f') 'N; ' num2str(hm0_struct.xyz(3),'%.1f') 'm' ')'];
        wshub_struct.xyz_str =  ['(' num2str(-1*hm0_struct.xyz(1),'%.3f') 'W; ' num2str(hm0_struct.xyz(2),'%.3f') 'N; ' num2str(hm0_struct.xyz(3),'%.1f') 'm' ')'];

        % do statistics
        [fitvals] = m_scatter(wshub_struct,hm0_struct,'density','quantiles',[0.05 0.5 0.95],'fit_func_Q','Poly','nolegend');
        save('fitvals.mat','fitvals');
        cd('../../');

    end

end

%% Hm0/Tp Directional

params1 = {'Hm0_Total','Hm0_Sea','Hm0_Swell'};
params2 = {'WindSpd_Hub'};
direc1 = {'MWD_Total','MWD_Sea','MWD_Swell'};
bins_hm0 = {0:2:8, 0:2:8, 0:2:8};
bins_wshub = {0:5:45};

cd([odir_hmo_wshub 'Directional']);

for i=1:no_locs

    % load structure
    fname = [fdir 'ASOW' num2str(i) '_all_structs.mat'];
    load(fname);

    for p = 1:length(params1)

        mkdir(params1{p});
        cd(params1{p});
        mkdir(['ASOW' num2str(i)]);
        cd(['ASOW' num2str(i)])
        
        hm0_struct = asow_params.(params1{p});
        hm0_struct.data(hm0_struct.data<0)=NaN;
        wshub_struct = asow_params.(params2{1});
        wshub_struct.data(wshub_struct.data<0)=NaN;

        mwd_struct = asow_params.(direc1{p});
        mwd_struct.data(mwd_struct.data<0)=NaN;

        % modify struct items for plotting
        hm0_struct.bins = bins_hm0{p};
        wshub_struct.bins = bins_wshub{1};
        hm0_struct.legend = 'SW_{US-EC}';
        wshub_struct.legend = 'SW_{US-EC}';
        wshub_struct.label = 'WS_{Hub}';
        wshub_struct.unit = 'm/s';
        hm0_struct.ttt_str_long = ' (1979-01-15–2021-12-31; \Deltat=1h) ';
        hm0_struct.ttt = [datenum('1979-01-15') datenum('2021-12-31') 60];
        hm0_struct.xyz_str =  ['(' num2str(-1*hm0_struct.xyz(1),'%.3f') 'W; ' num2str(hm0_struct.xyz(2),'%.3f') 'N; ' num2str(hm0_struct.xyz(3),'%.1f') 'm' ')'];
        wshub_struct.xyz_str =  ['(' num2str(-1*hm0_struct.xyz(1),'%.3f') 'W; ' num2str(hm0_struct.xyz(2),'%.3f') 'N; ' num2str(hm0_struct.xyz(3),'%.1f') 'm' ')'];

        % do statistics
        [fitvals] = m_scatter(wshub_struct,hm0_struct,'directional',mwd_struct,'density','quantiles',[0.05 0.5 0.95],'fit_func_Q','Poly','nolegend');
        save('fitvals.mat','fitvals');
        cd('../../');

    end

end