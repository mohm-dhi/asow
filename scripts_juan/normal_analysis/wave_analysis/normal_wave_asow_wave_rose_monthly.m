%%
clear all; clc

addpath(genpath('C:\DHI\Projects\AtlanticShores\Scripts\asow\scripts_juan\normal_analysis'));


%%

fdir = 'C:\DHI\Projects\AtlanticShores\Data\TimeSeries\Deliverable\CSVFinal\structs\';
odir_wr = 'C:\DHI\Projects\AtlanticShores\Analysis\Normal\Waves\Results\WaveRose\';
main_dir = 'Hm0';
no_locs = 7;


%% Wave yeights

params = {'Hm0_Total','Hm0_Sea','Hm0_Swell'};
direc = {'MWD_Total','MWD_Sea','MWD_Swell'};
bins = {0:3:12, 0:2:10, 0:2:10};
leg = {'H_{m0,Total}','H_{m0,Sea}','H_{m0,Swell}'};
di = {[0 0.01 0.1 0.25 0.5 0.75 1:0.5:4],[0 0.01 0.1 0.25:0.25:1 1:0.5:4],[0 0.05 0.1 0.25:0.25:1 1:0.25:2 2:0.5:4]};

for i=1:no_locs

    cd([odir_wr 'P' num2str(i)]);
    mkdir(main_dir);
    mkdir(main_dir,'Monthly');
    cd([main_dir '\Monthly']);

    % load structure
    fname = [fdir 'ASOW' num2str(i) '_all_structs.mat'];
    load(fname);

    for p = 1:length(params)

        mkdir(params{p});
        cd(params{p});

        wave_struct = asow_params.(params{p});
        wave_struct.data(wave_struct.data<0)=NaN;
        mwd_struct = asow_params.(direc{p});
        mwd_struct.data(mwd_struct.data<0)=NaN;

        % modify struct items for plotting
        wave_struct.bins = bins{p};
        %     wave_struct.unit = 'mMLLW';
        wave_struct.legend = 'SW_{US-EC}';
        wave_struct.ttt_str_long = ' (1979-01-15–2021-12-31; \Deltat=1h) ';
        wave_struct.ttt = [datenum('1979-01-15') datenum('2021-12-31') 60];
        wave_struct.xyz_str =  ['(' num2str(-1*wave_struct.xyz(1),'%.3f') 'W; ' num2str(wave_struct.xyz(2),'%.3f') 'N; ' num2str(wave_struct.xyz(3),'%.1f') 'm' ')'];

        % do statistics
        m_scatter(wave_struct,mwd_struct,'monthly','di',di{p});
        cd('..');

    end

end