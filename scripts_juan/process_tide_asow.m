%% load paths
addpath(genpath('C:\Users\jngz\MATLAB\Projects\potlab\src'));
addpath(genpath('C:\Users\jngz\MATLAB\Projects\potlab\res'));
% addpath('C:\Users\jngz\OneDrive - DHI\2022\41806529-AtlanticShores\Scripts\FileIO');

%%
fdir = 'C:\DHI\Projects\AtlanticShores\Data\TimeSeries\Points to be delivered SW_HD\';
fname = {'P1_US_EastCoast_HD_-73.945_39.651_24.6_3051.4_1979-01-01_2021-12-31_.csv';
    'P2_US_EastCoast_HD_-73.95_39.307_29.9_3458.2_1979-01-01_2021-12-31_.csv';
    'P3_US_EastCoast_HD_-74.044_39.201_24.4_2725.6_1979-01-01_2021-12-31_.csv';
    'P4_US_EastCoast_HD_-74.116_39.161_28.5_2678.6_1979-01-01_2021-12-31_.csv';
    'P5_US_EastCoast_HD_-74.118_39.252_26.8_3106.9_1979-01-01_2021-12-31_.csv';
    'P6_US_EastCoast_HD_-74.111_39.357_23.2_3174_1979-01-01_2021-12-31_.csv';
    'P7_US_EastCoast_HD_-74.212_39.284_21.9_3302_1979-01-01_2021-12-31_.csv'};
xyh = [-73.945 39.651 24.6;
    -73.95 39.307 29.9;
    -74.044 39.201 24.4;
    -74.116 39.161 28.5;
    -74.118 39.252 26.8;
    -74.111 39.357 23.2;
    -74.212 39.284 21.9];
headerlines = 16;
tmin            = datenum([1979 01 01 01 00 00]);  % Full HD period start
tmax            = datenum([2021 12 31 23 00 00]);  % Full HD period end
tttt_CURR       = [tmin tmax  30];
bins_WL         = [-2:0.25:2];
lgnd_hd_us  = 'HD_{US-EC}';

%% process water levels

for i=1:length(fname)

    name = ['ASOW' num2str(i)];
    mkdir(name);
    cd(name);

    % read water levels from MOOD
    [wlevel_arr,curr2d_arr] = read_mood_hd_csv(fdir,fname{i},headerlines);

    % WL from m_tide
    asow_wl = m_structure(name,xyh(i,:),tttt_CURR,lgnd_hd_us,wlevel_arr,['WL'],1,bins_WL);

    % run u-tide
    m_tide(asow_wl)

    cd('..');

end

%% process speed and direction

for i=1:length(fname)

    name = ['ASOW' num2str(i)];
    mkdir(name);
    cd(name);

    [wlevel_arr,curr2d_arr] = read_mood_hd_csv(fdir,fname{i},headerlines);

    asow_u    = m_structure(name,xyh(i,:),tttt_CURR,lgnd_hd_us,curr2d_arr,['uvel'],1,-1:0.01:1);
    asow_v    = m_structure(name,xyh(i,:),tttt_CURR,lgnd_hd_us,curr2d_arr,['vvel'],2,-1:0.01:1);
    [asow_u.data, asow_v.data] = spddir2uv(asow_u.data,asow_v.data); % uv2spddir
    
    asow_u.unit='m/s';
    asow_v.unit='m/s';

    m_tide(asow_u,asow_v);

    cd('..');

end