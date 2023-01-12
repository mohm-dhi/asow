%%
clear all; clc

addpath(genpath('C:\DHI\Projects\AtlanticShores\Scripts\asow\scripts_juan\normal_analysis'));


%%

fdir = 'C:\DHI\Projects\AtlanticShores\Data\TimeSeries\Deliverable\CSVFinal\structs\';
odir_m = 'C:\DHI\Projects\AtlanticShores\Analysis\Normal\Waves\DirectionalStats\';
no_locs = 7;


%% Wave yeights, monthly

params = {'Hm0_Total','Hm0_Sea','Hm0_Swell'};
direc = {'MWD_Total','MWD_Sea','MWD_Swell'};
bins = {0:3:12, 0:2:10, 0:2:10};
leg = {'Hm0_{Total}','Hm0_{Sea}','Hm0_{Swell}'};

for i=1:no_locs

    cd(odir_m);

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
        mwd_struct = asow_params.(direc{p});
        mwd_struct.data(mwd_struct.data<0)=NaN;

        % modify struct items for plotting
        wave_struct.bins = bins{p};
        %     wave_struct.unit = 'mMLLW';
        wave_struct.legend = leg{p};
        wave_struct.ttt_str_long = ' (1979-01-15–2021-12-31; \Deltat=1h) ';
        wave_struct.ttt = [datenum('1979-01-15') datenum('2021-12-31') 60];
        wave_struct.xyz_str =  ['(' num2str(-1*wave_struct.xyz(1),'%.3f') 'W;' num2str(wave_struct.xyz(2),'%.3f') 'N;' num2str(wave_struct.xyz(3),'%.1f') 'm' ')'];

        % do statistics
        m_statistics(wave_struct,'directional',mwd_struct);
        cd('../../');

    end

end


%% Wave Tp, monthly

params = {'Tp_Total','Tp_Sea','Tp_Swell'};
leg = {'T_{p,Total}','T_{p,Sea}','T_{p,Swell}'};
bins = {0:5:30, 0:5:25, 0:5:30};

for i=1:no_locs

    cd(odir_m);

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
        wave_struct.xyz_str =  ['(' num2str(-1*wave_struct.xyz(1),'%.3f') 'W;' num2str(wave_struct.xyz(2),'%.3f') 'N;' num2str(-1*wave_struct.xyz(3),'%.1f') 'm' ')'];

        % do statistics
        m_statistics(wave_struct,'monthly');
        cd('../../');

    end

end

%% Wave T01, monthly

params = {'T01_Total','T01_Sea','T01_Swell'};
leg = {'T_{01,Total}','T_{01,Sea}','T_{01,Swell}'};
bins = {0:5:25, 0:5:20, 0:5:25};

for i=1:no_locs

    cd(odir_m);

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
        wave_struct.xyz_str =  ['(' num2str(-1*wave_struct.xyz(1),'%.3f') 'W;' num2str(wave_struct.xyz(2),'%.3f') 'N;' num2str(-1*wave_struct.xyz(3),'%.1f') 'm' ')'];

        % do statistics
        m_statistics(wave_struct,'monthly');
        cd('../../');

    end

end

%% Wave T01, monthly

params = {'T02_Total','T02_Sea','T02_Swell'};
leg = {'T_{02,Total}','T_{02,Sea}','T_{02,Swell}'};
bins = {0:5:20, 0:5:20, 0:5:20};

for i=1:no_locs

    cd(odir_m);

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
        wave_struct.xyz_str =  ['(' num2str(-1*wave_struct.xyz(1),'%.3f') 'W;' num2str(wave_struct.xyz(2),'%.3f') 'N;' num2str(-1*wave_struct.xyz(3),'%.1f') 'm' ')'];

        % do statistics
        m_statistics(wave_struct,'monthly');
        cd('../../');

    end

end