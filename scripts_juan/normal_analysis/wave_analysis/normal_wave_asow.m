%%
addpath(genpath('C:\Users\jngz\MATLAB\Projects\potlab\src'));
addpath(genpath('C:\Users\jngz\MATLAB\Projects\potlab\res'));
addpath('C:\DHI\Projects\AtlanticShores\Scripts\asow\scripts_juan\normal_analysis\tidal_analysis');
addpath('C:\DHI\Projects\AtlanticShores\Scripts\asow\scripts_juan\normal_analysis');


%%

fdir = 'C:\DHI\Projects\AtlanticShores\Data\TimeSeries\Structs\';
odir_m = 'C:\DHI\Projects\AtlanticShores\Analysis\Normal\Waves\MonthlyStats\';
odir_y = 'C:\DHI\Projects\AtlanticShores\Analysis\Normal\Waves\AnnualStats\';
no_locs = 7;

mllw = [-0.55,-0.58,-0.60,-0.61,-0.61,-0.60,-0.63]; %vdatum

%% Wave yeights, yearly

params = {'Hm0_Total','Hm0_Sea','Hm0_Swell'};
leg = {'Hm0_{Total}','Hm0_{Sea}','Hm0_{Swell}'};
bins = {0:3:12, 0:2:10, 0:2:10};

for i=1:no_locs

    cd(odir_y);
    mkdir(['ASOW' num2str(i)]);
    cd(['ASOW' num2str(i)])

    % load structure
    fname = [fdir 'ASOW' num2str(i) '_all_structs.mat'];
    load(fname);

    for p = 1:length(params)

        % calculate tidal parameters
        wave_struct = asow_params_out.(params{p});

        % modify struct items for plotting
        wave_struct.bins = bins{p};
        %     wave_struct.unit = 'mMLLW';
        wave_struct.legend = leg{p};
        wave_struct.ttt_str_long = ' (1979-01-15–2021-12-31; \Deltat=60min) ';
        wave_struct.ttt = [datenum('1979-01-15') datenum('2021-12-31') 60];
        %     wave_struct.xyz_str =  ['('  num2str(-1*wave_struct.xyz(1)) 'W;' num2str(wave_struct.xyz(2)) 'N;'  num2str(-1*(wave_struct.xyz(3)+abs(mllw(i)))) wave_struct.unit ')'];

        % do statistics
        m_statistics(wave_struct,'yearly');

    end

end

%% Wave yeights, monthly

params = {'Hm0_Total','Hm0_Sea','Hm0_Swell'};
bins = {0:3:12, 0:2:10, 0:2:10};
leg = {'Hm0_{Total}','Hm0_{Sea}','Hm0_{Swell}'};

for i=1:no_locs

    cd(odir_m);
    mkdir(['ASOW' num2str(i)]);
    cd(['ASOW' num2str(i)])

    % load structure
    fname = [fdir 'ASOW' num2str(i) '_all_structs.mat'];
    load(fname);

    for p = 1:length(params)

        % calculate tidal parameters
        wave_struct = asow_params_out.(params{p});

        % modify struct items for plotting
        wave_struct.bins = bins{p};
        %     wave_struct.unit = 'mMLLW';
        wave_struct.legend = leg{p};
        wave_struct.ttt_str_long = ' (1979-01-15–2021-12-31; \Deltat=60min) ';
        wave_struct.ttt = [datenum('1979-01-15') datenum('2021-12-31') 60];
        %     wave_struct.xyz_str =  ['('  num2str(-1*wave_struct.xyz(1)) 'W;' num2str(wave_struct.xyz(2)) 'N;'  num2str(-1*(wave_struct.xyz(3)+abs(mllw(i)))) wave_struct.unit ')'];

        % do statistics
        m_statistics(wave_struct,'monthly');

    end

end