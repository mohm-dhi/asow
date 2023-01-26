%%
clear all; clc;

% addpath(genpath('C:\Users\jngz\MATLAB\Projects\potlab\res'));
% addpath(genpath('C:\Users\jngz\MATLAB\Projects\potlab\src'));
% addpath('C:\DHI\Projects\AtlanticShores\Scripts\asow\scripts_juan\normal_analysis\tidal_analysis');
addpath(genpath('C:\DHI\Projects\AtlanticShores\Scripts\asow\scripts_juan\normal_analysis'));

%%

fdir = 'C:\DHI\Projects\AtlanticShores\Data\TimeSeries\TidalAnalysis\wl\MSL\';
odir_t = 'C:\DHI\Projects\AtlanticShores\Analysis\Normal\Tidal\TidalLevels\';
odir_y = 'C:\DHI\Projects\AtlanticShores\Analysis\Normal\Tidal\AnnualStats\';
odir_o = 'C:\DHI\Projects\AtlanticShores\Analysis\Normal\Tidal\Omni\';
no_locs = 7;

mllw = [-0.55,-0.58,-0.60,-0.61,-0.61,-0.60,-0.63]; %vdatum

%% Tidal Levels (wrt MSL)

cd(odir_t);

for i=1:no_locs

    % load structure
    fname = [fdir 'ASOW' num2str(i) '\ASOW' num2str(i) '_wl.mat'];
    load(fname);

    % find analysis time, apply to data
    ti = find(WLTde.datetime>=datetime(1979,1,15,0,0,0));
    WLTde.time = WLTde.time(ti);
    WLTde.datetime = WLTde.datetime(ti);
    WLTde.data = WLTde.data(ti);

    % calculate annual stats
%     WL_Tide = asow_params.WL_Tide;
%     WL_Tide.data = WL_Tide.data-abs(mllw(i));
    tidal_levels = m_tidal_levels(WLTde);

    % min and max 
    tidal_levels.names_c = {'min','max'};
    tidal_levels.values_c = [min(WLTde.data) max(WLTde.data)];

    % save to file
    save(['MSL\asow' num2str(i) '_tidal_levels_MSL.mat'], 'tidal_levels');

end

%% Tidal Levels (wrt MLLW)

% cd(odir_t);
% 
% for i=1:no_locs
% 
%     % load structure
%     fname = [fdir 'ASOW' num2str(i) '\ASOW' num2str(i) '_wl.mat'];
%     load(fname);
% 
%     % find analysis time, apply to data
%     ti = find(WLTde.datetime>=datetime(1979,1,15,0,0,0));
%     WLTde.time = WLTde.time(ti);
%     WLTde.datetime = WLTde.datetime(ti);
%     WLTde.data = WLTde.data(ti)+abs(mllw(i));
% 
%     % calculate tidal parameters
% %     WL_Tide = asow_params.WL_Tide;
%     tidal_levels = m_tidal_levels(WLTde);
% 
%     % min and max 
%     tidal_levels.names_c = {'min','max'};
%     tidal_levels.values_c = [min(WLTde.data) max(WLTde.data)];
% 
%     % save to file
%     save(['MLLW\asow' num2str(i) '_tidal_levels_MLLW.mat'], 'tidal_levels');
% 
% end

%% Calculate yearly statistics on Total WL (wrt MSL) 

for i=1:no_locs

    cd([odir_y 'MSL']);
    mkdir(['P' num2str(i)]);
    cd(['P' num2str(i)])

    % load structure
    fname = [fdir 'ASOW' num2str(i) '\ASOW' num2str(i) '_wl.mat'];
    load(fname);

    % find analysis time, apply to data
    ti = find(WL.datetime>=datetime(1979,1,15,0,0,0));
    WL.time = WL.time(ti);
    WL.datetime = WL.datetime(ti);
    WL.data = WL.data(ti);
    WL.name = ['P' num2str(i)]

    % modify struct items for plotting
    WL.bins = -2:0.5:3;
    WL.unit = 'mMSL';
    WL.legend = 'WLTotal MSL';
    WL.ttt_str_long = ' (1979-01-15–2021-12-31; \Deltat=30min) ';
    WL.ttt = [datenum('1979-01-15') datenum('2021-12-31') 60];
    WL.xyz_str =  [' ('  num2str(-1*WL.xyz(1),'%.3f') 'W; ' num2str(WL.xyz(2),'%.3f') 'N; '  num2str(-1*WL.xyz(3),'%.1f') WL.unit ') '];

    % do statistics
    m_statistics(WL,'yearly');

end

%% Calculate yearly statistics on Total WL (wrt MLLW) 

% for i=1:no_locs
% 
%     cd([odir_y 'MLLW']);
%     mkdir(['ASOW' num2str(i)]);
%     cd(['ASOW' num2str(i)])
% 
%     % load structure
%     fname = [fdir 'ASOW' num2str(i) '\ASOW' num2str(i) '_wl.mat'];
%     load(fname);
% 
%     % find analysis time, apply to data
%     ti = find(WL.datetime>=datetime(1979,1,15,0,0,0));
%     WL.time = WL.time(ti);
%     WL.datetime = WL.datetime(ti);
%     WL.data = WL.data(ti)+abs(mllw(i));
% 
% 
%     % modify struct items for plotting
%     WL.bins = -1.2:0.5:4;
%     WL.unit = 'mMSL';
%     WL.legend = 'WLTotal MLLW';
%     WL.ttt_str_long = ' (1979-01-15–2021-12-31; \Deltat=30min) ';
%     WL.ttt = [datenum('1979-01-15') datenum('2021-12-31') 60];
%     WL.xyz_str =  ['('  num2str(-1*WL.xyz(1)) 'W;' num2str(WL.xyz(2)) 'N;'  num2str(-1*(WL.xyz(3))) WL.unit ')'];
% 
%     % do statistics
%     m_statistics(WL,'yearly');
% 
% end