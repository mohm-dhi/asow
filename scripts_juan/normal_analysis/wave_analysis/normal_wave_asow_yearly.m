%%
clear all; clc

addpath(genpath('C:\DHI\Projects\AtlanticShores\Scripts\asow\scripts_juan\normal_analysis'));


%%

fdir = 'C:\DHI\Projects\AtlanticShores\Data\TimeSeries\Deliverable\CSVFinal\structs\';
odir_y = 'C:\DHI\Projects\AtlanticShores\Analysis\Normal\Waves\AnnualStats\';
no_locs = 7;


%% Wave yeights, yearly

params = {'Hm0_Total','Hm0_Sea','Hm0_Swell'};
bins = {0:3:12, 0:2:10, 0:2:10};
leg = {'Hm0_{Total}','Hm0_{Sea}','Hm0_{Swell}'};

for i=1:no_locs

    cd(odir_y);

    % load structure
    fname = [fdir 'ASOW' num2str(i) '_all_structs.mat'];
    load(fname);

    for p = 1:length(params)

        mkdir(params{p});
        cd(params{p});
        mkdir(['ASOW' num2str(i)]);
        cd(['ASOW' num2str(i)])
        wave_struct = asow_params.(params{p});
        wave_struct.data(wave_struct.data<0)=NaN;

        % modify struct items for plotting
        wave_struct.bins = bins{p};
        %     wave_struct.unit = 'mMLLW';
        wave_struct.legend = leg{p};
        wave_struct.ttt_str_long = ' (1979-01-15–2021-12-31; \Deltat=1h) ';
        wave_struct.ttt = [datenum('1979-01-15') datenum('2021-12-31') 60];
        wave_struct.xyz_str =  ['(' num2str(-1*wave_struct.xyz(1),'%.3f') 'W;' num2str(wave_struct.xyz(2),'%.3f') 'N;' num2str(wave_struct.xyz(3),'%.1f') 'm' ')'];

        % do statistics
        m_statistics(wave_struct,'yearly');
        cd('../../');

    end

end


%% Wave Tp, yearly

params = {'Tp_Total','Tp_Sea','Tp_Swell'};
leg = {'T_{p,Total}','T_{p,Sea}','T_{p,Swell}'};
bins = {0:5:30, 0:5:25, 0:5:30};

for i=1:no_locs

    cd(odir_y);

    % load structure
    fname = [fdir 'ASOW' num2str(i) '_all_structs.mat'];
    load(fname);

    for p = 1:length(params)

        mkdir(params{p});
        cd(params{p});
        mkdir(['ASOW' num2str(i)]);
        cd(['ASOW' num2str(i)])
        wave_struct = asow_params.(params{p});
        wave_struct.data(wave_struct.data<0)=NaN;

        % modify struct items for plotting
        wave_struct.bins = bins{p};
        %     wave_struct.unit = 'mMLLW';
        wave_struct.legend = leg{p};
        wave_struct.ttt_str_long = ' (1979-01-15–2021-12-31; \Deltat=1h) ';
        wave_struct.ttt = [datenum('1979-01-15') datenum('2021-12-31') 60];
        wave_struct.xyz_str =  ['(' num2str(-1*wave_struct.xyz(1),'%.3f') 'W;' num2str(wave_struct.xyz(2),'%.3f') 'N;' num2str(wave_struct.xyz(3),'%.1f') 'm' ')'];

        % do statistics
        m_statistics(wave_struct,'yearly');
        cd('../../');

    end

end

%% Wave T01, yearly

params = {'T01_Total','T01_Sea','T01_Swell'};
leg = {'T_{01,Total}','T_{01,Sea}','T_{01,Swell}'};
bins = {0:5:25, 0:5:20, 0:5:25};

for i=1:no_locs

    cd(odir_y);

    % load structure
    fname = [fdir 'ASOW' num2str(i) '_all_structs.mat'];
    load(fname);

    for p = 1:length(params)

        mkdir(params{p});
        cd(params{p});
        mkdir(['ASOW' num2str(i)]);
        cd(['ASOW' num2str(i)])
        wave_struct = asow_params.(params{p});
        wave_struct.data(wave_struct.data<0)=NaN;

        % modify struct items for plotting
        wave_struct.bins = bins{p};
        %     wave_struct.unit = 'mMLLW';
        wave_struct.legend = leg{p};
        wave_struct.ttt_str_long = ' (1979-01-15–2021-12-31; \Deltat=1h) ';
        wave_struct.ttt = [datenum('1979-01-15') datenum('2021-12-31') 60];
        wave_struct.xyz_str =  ['(' num2str(-1*wave_struct.xyz(1),'%.3f') 'W;' num2str(wave_struct.xyz(2),'%.3f') 'N;' num2str(wave_struct.xyz(3),'%.1f') 'm' ')'];

        % do statistics
        m_statistics(wave_struct,'yearly');
        cd('../../');

    end

end

%% Wave T01, yearly

params = {'T02_Total','T02_Sea','T02_Swell'};
leg = {'T_{02,Total}','T_{02,Sea}','T_{02,Swell}'};
bins = {0:5:20, 0:5:20, 0:5:20};

for i=1:no_locs

    cd(odir_y);

    % load structure
    fname = [fdir 'ASOW' num2str(i) '_all_structs.mat'];
    load(fname);

    for p = 1:length(params)

        mkdir(params{p});
        cd(params{p});
        mkdir(['ASOW' num2str(i)]);
        cd(['ASOW' num2str(i)])
        wave_struct = asow_params.(params{p});
        wave_struct.data(wave_struct.data<0)=NaN;

        % modify struct items for plotting
        wave_struct.bins = bins{p};
        %     wave_struct.unit = 'mMLLW';
        wave_struct.legend = leg{p};
        wave_struct.ttt_str_long = ' (1979-01-15–2021-12-31; \Deltat=1h) ';
        wave_struct.ttt = [datenum('1979-01-15') datenum('2021-12-31') 60];
        wave_struct.xyz_str =  ['(' num2str(-1*wave_struct.xyz(1),'%.3f') 'W;' num2str(wave_struct.xyz(2),'%.3f') 'N;' num2str(wave_struct.xyz(3),'%.1f') 'm' ')'];

        % do statistics
        m_statistics(wave_struct,'yearly');
        cd('../../');

    end

end