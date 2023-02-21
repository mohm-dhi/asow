%%
clear all; clc

addpath(genpath('C:\DHI\Projects\AtlanticShores\Scripts\asow\scripts_juan\normal_analysis'));

%%

fdir = 'C:\DHI\Projects\AtlanticShores\Data\TimeSeries\Deliverable\CSVFinal\structs\';
odir_scatter = 'C:\DHI\Projects\AtlanticShores\Analysis\Normal\Waves\Results\Scatter\';
main_dir = 'WDir_MWD';
no_locs = 7;

%% WindDir/MWD Omni

params1 = {'WindDir_10m'};
params2 = {'MWD_Total','MWD_Sea','MWD_Swell'};
bins_winddir = {0:30:330};
bins_mwd = {0:30:330};

for i=1:no_locs

    cd([odir_scatter 'P' num2str(i)]);
    mkdir(main_dir);
    mkdir(main_dir,'Omni');
    cd([main_dir '\Omni']);

    % load structure
    fname = [fdir 'ASOW' num2str(i) '_all_structs.mat'];
    load(fname);

    for p = 1:length(params2)

        mkdir(params2{p});
        cd(params2{p});

        winddir_struct = asow_params.(params1{1});
        winddir_struct.data(winddir_struct.data<0)=NaN;
        mwd_struct = asow_params.(params2{p});
        mwd_struct.data(mwd_struct.data<0)=NaN;

        % modify struct items for plotting
        winddir_struct.bins = bins_winddir{1};
        mwd_struct.bins = bins_mwd{1};
        winddir_struct.legend = 'SW_{US-EC}';
        winddir_struct.label = 'Wind Dir_{10m}';
        winddir_struct.unit = '\circN-from';
        mwd_struct.legend = 'SW_{US-EC}';
        winddir_struct.ttt_str_long = ' (1979-01-15–2021-12-31; \Deltat=1h) ';
        winddir_struct.ttt = [datenum('1979-01-15') datenum('2021-12-31') 60];
        winddir_struct.xyz_str =  ['(' num2str(-1*winddir_struct.xyz(1),'%.3f') 'W; ' num2str(winddir_struct.xyz(2),'%.3f') 'N; ' num2str(winddir_struct.xyz(3),'%.1f') 'm' ')'];

        % do scatter
        m_scatter(winddir_struct,mwd_struct,'monthly','plottype','scatter','density','nolegend');
        cd('..');

    end

end