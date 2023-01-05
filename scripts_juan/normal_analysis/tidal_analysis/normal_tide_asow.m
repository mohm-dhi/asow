%%
addpath(genpath('C:\Users\jngz\MATLAB\Projects\potlab\src'));
addpath(genpath('C:\Users\jngz\MATLAB\Projects\potlab\res'));
addpath('C:\DHI\Projects\AtlanticShores\Scripts\asow\scripts_juan\normal_analysis\tidal_analysis');
addpath('C:\DHI\Projects\AtlanticShores\Scripts\asow\scripts_juan\normal_analysis');


%%

fdir = 'C:\DHI\Projects\AtlanticShores\Data\TimeSeries\Structs\';
odir_t = 'C:\DHI\Projects\AtlanticShores\Analysis\Normal\Tidal\TidalLevels\';
odir_y = 'C:\DHI\Projects\AtlanticShores\Analysis\Normal\Tidal\AnnualStats\';
no_locs = 7;

mllw = [-0.55,-0.58,-0.60,-0.61,-0.61,-0.60,-0.63]; %vdatum

%% Tidal Levels (wrt MSL, ie remove MLLW correction)

cd(odir_t);

for i=1:no_locs

    % load structure
    fname = [fdir 'ASOW' num2str(i) '_all_structs.mat'];
    load(fname);

    % calculate annual stats
    WL_Tide = asow_params_out.WL_Tide;
    WL_Tide.data = WL_Tide.data-abs(mllw(i));
    tidal_levels = m_tidal_levels(WL_Tide);

    % min and max 
    tidal_levels.names_c = {'min','max'};
    tidal_levels.values_c = [min(WL_Tide.data) max(WL_Tide.data)];

    % save to file
    save(['asow' num2str(i) '_tidal_levels_MSL.mat'], 'tidal_levels');

end

%% Tidal Levels (wrt MLLW)

cd(odir_t);

for i=1:no_locs

    % load structure
    fname = [fdir 'ASOW' num2str(i) '_all_structs.mat'];
    load(fname);

    % calculate tidal parameters
    WL_Tide = asow_params_out.WL_Tide;
    tidal_levels = m_tidal_levels(WL_Tide);

    % min and max 
    tidal_levels.names_c = {'min','max'};
    tidal_levels.values_c = [min(WL_Tide.data) max(WL_Tide.data)];

    % save to file
    save(['asow' num2str(i) '_tidal_levels_MLLW.mat'], 'tidal_levels');

end

%% Calculate yearly statistics on Total WL (wrt MSL) 

for i=1:no_locs

    cd([odir_y 'MSL']);
    mkdir(['ASOW' num2str(i)]);
    cd(['ASOW' num2str(i)])

    % load structure
    fname = [fdir 'ASOW' num2str(i) '_all_structs.mat'];
    load(fname);

    % calculate tidal parameters
    WL_Total = asow_params_out.WL_Total;
    WL_Total.data = WL_Total.data-abs(mllw(i));

    % modify struct items for plotting
    WL_Total.bins = -2:0.5:3;
    WL_Total.unit = 'mMSL';
    WL_Total.legend = 'WLTotal MSL';
    WL_Total.ttt_str_long = ' (1979-01-15–2021-12-31; \Deltat=60min) ';
    WL_Total.ttt = [datenum('1979-01-15') datenum('2021-12-31') 60];
    WL_Total.xyz_str =  ['('  num2str(-1*WL_Total.xyz(1)) 'W;' num2str(WL_Total.xyz(2)) 'N;'  num2str(-1*WL_Total.xyz(3)) WL_Total.unit ')'];

    % do statistics
    m_statistics(WL_Total,'yearly');

end

%% Calculate yearly statistics on Total WL (wrt MLLW) 

for i=1:no_locs

    cd([odir_y 'MLLW']);
    mkdir(['ASOW' num2str(i)]);
    cd(['ASOW' num2str(i)])

    % load structure
    fname = [fdir 'ASOW' num2str(i) '_all_structs.mat'];
    load(fname);

    % calculate tidal parameters
    WL_Total = asow_params_out.WL_Total;

    % modify struct items for plotting
    WL_Total.bins = -1.2:0.5:4;
    WL_Total.unit = 'mMLLW';
    WL_Total.legend = 'WLTotal MLLW';
    WL_Total.ttt_str_long = ' (1979-01-15–2021-12-31; \Deltat=60min) ';
    WL_Total.ttt = [datenum('1979-01-15') datenum('2021-12-31') 60];
    WL_Total.xyz_str =  ['('  num2str(-1*WL_Total.xyz(1)) 'W;' num2str(WL_Total.xyz(2)) 'N;'  num2str(-1*(WL_Total.xyz(3)+abs(mllw(i)))) WL_Total.unit ')'];

    % do statistics
    m_statistics(WL_Total,'yearly');

end