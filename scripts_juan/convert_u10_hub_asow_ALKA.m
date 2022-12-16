addpath(genpath('\\USDEN1-STOR.DHI.DK\Projects\41806529\_metocean_scripts\potlab_v2\src'));
addpath(genpath('\\USDEN1-STOR.DHI.DK\Projects\41806529\_metocean_scripts\potlab_v2\res'));

% fdir_cfsr = 'C:\Users\jngz\OneDrive - DHI\2022\41806529-AtlanticShores\Data\CFSR MOOD\OutputLocations\';

fdir_cfsr = '\\USDEN1-STOR.DHI.DK\Projects\41806529\02_RAW_MOOD_data\CFSR interpolated\concatenated\';
outputdir = '\\USDEN1-STOR.DHI.DK\Projects\41806529\03_PROCESSED_MOOD_data\02-OutputLocations\CFSR_Wind\';

% cfsr_names = {'P1_Global_Wind_CFSR_286.055_39.651_10_22716.1_1979-01-01_2022-10-01_.dfs0';
%     'P2_Global_Wind_CFSR_286.05_39.307_10_22716.1_1979-01-01_2022-10-01_.dfs0';
%     'P3_Global_Wind_CFSR_285.956_39.201_10_22716.1_1979-01-01_2022-10-01_.dfs0';
%     'P4_Global_Wind_CFSR_285.884_39.161_10_22716.1_1979-01-01_2022-10-01_.dfs0';
%     'P5_Global_Wind_CFSR_285.882_39.252_10_22716.1_1979-01-01_2022-10-01_.dfs0';
%     'P6_Global_Wind_CFSR_285.889_39.357_10_22716.1_1979-01-01_2022-10-01_.dfs0';
%     'P7_Global_Wind_CFSR_285.788_39.284_10_22716.1_1979-01-01_2022-10-01_.dfs0'};

cfsr_names = {'P1_CFSR-global-_1979-01-01_-_2022-12-31.dfs0';
    'P2_CFSR-global-_1979-01-01_-_2022-12-31.dfs0';
    'P3_CFSR-global-_1979-01-01_-_2022-12-31.dfs0';
    'P4_CFSR-global-_1979-01-01_-_2022-12-31.dfs0';
    'P5_CFSR-global-_1979-01-01_-_2022-12-31.dfs0';
    'P6_CFSR-global-_1979-01-01_-_2022-12-31.dfs0';
    'P7_CFSR-global-_1979-01-01_-_2022-12-31.dfs0';
    'P2D_CFSR-global-_1979-01-01_-_2022-12-31.dfs0'};

% lon lat depth
xyh = [-73.945 39.651 24.6;
    -73.95 39.307 29.9;
    -74.044 39.201 24.4;
    -74.116 39.161 28.5;
    -74.118 39.252 26.8;
    -74.111 39.357 23.2;
    -74.212 39.284 21.9];

% bins and dates
bins_D = 0:30:360; % wind dir
bins_WS = 0:3:33;  % wind speed
bins_U = -30:0.1:30;
bins_V = -30:0.1:30;
ttt = [datenum([1979 01 01 00 00 00]) datenum([2021 12 31 23 00 00]) 60];

cfsr_height = 10;  % cfsr wind height
hub_height =  147; % mMSL

% for each cfsr file
for i = 1:length(cfsr_names)-1

    % x y z
    xyz = xyh(i,:);

    % station name
    name = ['ASOW' num2str(i)];

    % cfsr name
    cfsr_path = [fdir_cfsr cfsr_names{i}];

    % load cfsr .dfs0
    % change to load u,v

    CFSR_U = m_structure(name , xyz , ttt , 'Data' , cfsr_path, 'u-wind' , 1, bins_U);
    CFSR_V = m_structure(name , xyz , ttt , 'Data' , cfsr_path, 'v-wind' , 2, bins_V);

    [WindSpeed,WindDir] = uv2spddir(CFSR_U, CFSR_V, 'from');

    WindSpd_10m = m_structure(name , xyz , ttt , 'Data' , [CFSR_U.time WindSpeed.data], 'WS-10m' , 1, bins_WS);
    WindDir_10m = m_structure(name , xyz , ttt , 'Data' , [CFSR_U.time WindDir.data], 'WD-10m' , 1, bins_D);
    WindDir_Hub = m_structure(name , xyz , ttt , 'Data' , [CFSR_U.time WindDir.data], 'WD-Hub' , 1, bins_D);

    % m_structure for hub speed
    WindSpd_Hub = m_structure(name , xyz , ttt , 'Data' , [CFSR_U.time WindSpeed.data], 'WS-Hub' , 1, bins_WS);

    % convert to hub height
    Uz = U102Uz(WindSpd_10m.data,hub_height,'power',0.111);
    WindSpd_Hub.data = Uz;

    % save struct
    save([outputdir name '_WindSpd_Hub.mat'],'WindSpd_Hub');
    save([outputdir name '_WindDir_Hub.mat'],'WindDir_10m');
    save([outputdir name '_WindSpd_10m.mat'],'WindSpd_10m');
    save([outputdir name '_WindDir_10m.mat'],'WindDir_10m');

end