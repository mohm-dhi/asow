%%
addpath(genpath('\\USDEN1-STOR.DHI.DK\Projects\41806529\_metocean_scripts\potlab_v2\src'));
%addpath(genpath('\\USDEN1-STOR.DHI.DK\Projects\41806529\_metocean_scripts\potlab_v2\res'));
addpath('\\usden1-stor.dhi.dk\projects\41806529\_metocean_scripts\potlab_v2\res\UTide\')
%%
fdir = '\\USDEN1-STOR.DHI.DK\Projects\41806529\02_RAW_MOOD_data\Points to be delivered SW_HD\';
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
tttt_CURR       = [tmin tmax  60  60];
bins_WL         = [-1.80:0.25:1.80];
lgnd_HDOS_mike2d  = 'HD_{MIKE2D}';

%% process water levels
if 0
    for i=1:length(fname)
        
        name = ['ASOW' num2str(i)];
        mkdir(name);
        cd(name);
        
        % read water levels from MOOD
        [wlevel_arr,curr2d_arr] = read_mood_hd_csv(fdir,fname{i},headerlines);
        
        % WL from m_tide
        asow_wl = m_structure(name,xyh(i,:),tttt_CURR,lgnd_HDOS_mike2d,wlevel_arr,['WL'],1,bins_WL);
        
        % run u-tide
        m_tide(asow_wl)
        
        cd('..');
        
    end
end
%% process speed and direction

for i=1:length(fname)
    
    name = ['ASOW' num2str(i)];
    
    [wlevel_arr,curr2d_arr] = read_mood_hd_csv(fdir,fname{i},headerlines);
    
    CS    = m_structure(name,xyh(i,:),tttt_CURR,lgnd_HDOS_mike2d,curr2d_arr,'CS',1,0:0.05:2);
    CD    = m_structure(name,xyh(i,:),tttt_CURR,lgnd_HDOS_mike2d,curr2d_arr,'CD',2,0:22.5:360);
    
    figure;
    plot(CS.data); hold on
    plot(CD.data);
    
    mkdir(name);
    cd(name);
    
    m_scatter(CD,CS);
    
    
    
    m_tide(CS,CD);
    
    figure;
    plot(CS.data); hold on
    plot(CSTde.data)
    plot(CSRes.data)
    
    
    %% CHECK
    fs = 2;
    x = CSRes.data;
    
    y = fft(x);
    
    n = length(x);          % number of samples
    f = (0:n-1)*(fs/n);     % frequency range
    power = abs(y).^2/n;    % power of the DFT
    
    figure
    
    subplot(211)
    semilogx(f,power)
    xlabel('Frequency')
    ylabel('Power')
    title('Residual')
    hold on
    
    hold on
    plot([1/12.75, 1/12.75],[0 80])
    
    fs = 2;
    x = CSTde.data;
    
    y = fft(x);
    
    n = length(x);          % number of samples
    f = (0:n-1)*(fs/n);     % frequency range
    power = abs(y).^2/n;    % power of the DFT
    
    subplot(212)
    semilogx(f,power)
    xlabel('Frequency')
    ylabel('Power')
    title('Tidal Currents')
    hold on
    plot([1/12.75, 1/12.75],[0 80])
    
    
    %ylim([0 1])
    
    cd('..');
    
end