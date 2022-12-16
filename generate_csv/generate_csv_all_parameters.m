% generate deliverable for Atlantic Shores
% Wave spectral data: ASOW_P1-DHI-MET-yyyy Metocean data_LocXXXX_x.0.csv
% Where CONTRACTOR is the Contractor acronym, yyyy the year of delivery
% XXXX is the location name and
% x.0 is the revision status.
%
% ALKA, Dec 2022

clear all; close all; clc

% Directories are:
% CFSR_Wind: contains the 10m (speed and direction) and hub height (speed) CFSR winds (one per location, per wind parameter, so three files per location)
% 
% TidalAnalysis\WL\MSL: contains the UTide results at MSL, each location has its own directory
% 
% TidalAnalysis\WL\MLLW: contains the total water level, tidal water level, and residual tidal level at MLLW, each location has its own directory
% 
% TidalAnalysis\UV: contains the speed and current tidal analysis for each location as a directory. Inside each directory there are four subdirs level1,...,level4 where level1 is bottom depth (100%), level2 is 75% depth, level3 is 50% depth and level4 is 25% depth


%% SETTINGS
% time series vector (to be uniform for all parameters)
time = (datetime(1979,1,15):hours(1):datetime(2021,12,31,23,0,0))';

% output directory
output_dir = '\\USDEN1-STOR.DHI.DK\Projects\41806529\07_Timeseries_CSV_Deliverables\';

% template
template_fname = '\\USDEN1-STOR.DHI.DK\Projects\41806529\07_Timeseries_CSV_Deliverables\TEMPLATE_TimeSeries';

% list of parameters
params_fname = '\\USDEN1-STOR.DHI.DK\Projects\41806529\07_Timeseries_CSV_Deliverables\ListParameters.xlsx';

% wave parameters directly from MOOD
SW_path = '\\USDEN1-STOR.DHI.DK\Projects\41806529\02_RAW_MOOD_data\Points to be delivered SW_HD\';

% wind paramters
wind_path = '\\USDEN1-STOR.DHI.DK\Projects\41806529\03_PROCESSED_MOOD_data\02-OutputLocations\CFSR_Wind\';

% water levels (post-processed)
wl_path = '\\USDEN1-STOR.DHI.DK\Projects\41806529\03_PROCESSED_MOOD_data\02-OutputLocations\TidalAnalysis\WL\MLLW\';

% currents (post-processed)
cur_path = '\\USDEN1-STOR.DHI.DK\Projects\41806529\03_PROCESSED_MOOD_data\02-OutputLocations\TidalAnalysis\UV\';

% current levels (following the order in the current path)
cur_z = {'Surf','75','Mid','25','Bot'};

%% READ TEMPLATE
template = readlines(template_fname);

% read overall parameters table
params = readtable(params_fname,'Sheet','Params');

% read how to translate MOOD header to final header
mood2tab = readtable(params_fname,'Sheet','TranslateWaves');

%% LOOP THROUGH OUTPUT LOCATIONS
for loc = 1:7
    
    %% INITIALIZE OUTPUT TABLE
    tab = array2table(nan(length(time),height(params)),...
        'VariableNames',params.Parameter);
    tab.Time = time;
    
    %% WAVES
    % find file
    sw_list = dir([SW_path,sprintf('P%d',loc),'*SW*']);
    
    % load file
    wave = readtable([sw_list.folder,'/',sw_list.name]);
    
    % ensure time matches and add to final table
    [~,Ia,Ib] = intersect(time,wave.datetime_ISO8601__UTC_);
    
    % translate MOOD header to final table
    for n = 1:height(mood2tab)
        tab.(mood2tab.Final{n})(Ia) = wave.(mood2tab.MOOD{n})(Ib);
    end
    
    %% WIND - 10m and hub height
    % load files
    load([wind_path,'ASOW',num2str(loc),'_CFSR_Spd.mat']);
    load([wind_path,'ASOW',num2str(loc),'_CFSR_Spd_Dir.mat']);
    load([wind_path,'ASOW',num2str(loc),'_CFSR_Spd_Hub.mat']);
    
    % ensure time matches and add to final table
    [~,Ia,Ib] = intersect(time,CFSR_Spd.datetime);
    
    tab.WindSpd_10m(Ia) = CFSR_Spd.data(Ib);
    tab.WindDir_10m(Ia) = CFSR_Dir.data(Ib);
    tab.WindSpd_Hub(Ia) = CFSR_Spd_Hub.data(Ib);
    tab.WindDir_Hub(Ia) = CFSR_Dir.data(Ib);
    
    %% WATER LEVELS
    % load files
    load([wl_path,'ASOW',num2str(loc),'/ASOW',num2str(loc),'_wl_mllw.mat']);
    load([wl_path,'ASOW',num2str(loc),'/ASOW',num2str(loc),'_wltide_mllw.mat']);
    load([wl_path,'ASOW',num2str(loc),'/ASOW',num2str(loc),'_wlres_mllw.mat']);
    
    % Total WL
    % ensure time matches and add to final table
    [~,Ia,Ib] = intersect(time,WL.datetime);
    tab.WL_Total(Ia) = WL.data(Ib);
 
    % Tidal WL
    % ensure time matches and add to final table
    [~,Ia,Ib] = intersect(time,WLTide.datetime);   
    tab.WL_Tide(Ia) = WLTide.data(Ib);
    
    % Residual WL
    % ensure time matches and add to final table
    [~,Ia,Ib] = intersect(time,WLRes.datetime);   
    tab.WL_Res(Ia) = WLRes.data(Ib);
    
    %% CURRENTS DEPTH AVERAGED
    load([cur_path,'ASOW',num2str(loc),'/ASOW',num2str(loc),'_current_analysis.mat']);
    % ensure time matches and add to final table
    [~,Ia,Ib] = intersect(time,CS.datetime);
    tab.CurSpd_DepAvg_Total(Ia) = CS.data(Ib);
    tab.CurDir_DepAvg_Total(Ia) = CD.data(Ib);
    
    tab.CurSpd_DepAvg_Tide(Ia) = CSTde.data(Ib);
    tab.CurDir_DepAvg_Tide(Ia) = CDTde.data(Ib);
    
    tab.CurSpd_DepAvg_Res(Ia) = CSRes.data(Ib);
    tab.CurDir_DepAvg_Res(Ia) = CDRes.data(Ib);   
    
    %% CURRENTS (PER DEPTH)
    for z = 1:length(cur_z)
        % load tiles
        % tidal and residual
        load([cur_path,'ASOW',num2str(loc),'/level_',num2str(z),'/U_total']);
        load([cur_path,'ASOW',num2str(loc),'/level_',num2str(z),'/D_total']);
        load([cur_path,'ASOW',num2str(loc),'/level_',num2str(z),'/U_res']);
        load([cur_path,'ASOW',num2str(loc),'/level_',num2str(z),'/D_res']);
        load([cur_path,'ASOW',num2str(loc),'/level_',num2str(z),'/U_tide']);
        load([cur_path,'ASOW',num2str(loc),'/level_',num2str(z),'/D_tide']);
        
        % create output varname depending on depth
        vdir = ['CurDir_',cur_z{z},'_Total'];
        vtotal = ['CurSpd_',cur_z{z},'_Total'];
        
        vdir_tidal = ['CurDir_',cur_z{z},'_Tide'];
        vtidal = ['CurSpd_',cur_z{z},'_Tide'];
        
        vdir_res = ['CurDir_',cur_z{z},'_Res'];
        vres = ['CurSpd_',cur_z{z},'_Res'];

        % Total Currents
        [~,Ia,Ib] = intersect(time,U_total.datetime);
        tab.(vdir)(Ia) = D_total.data(Ib);
        tab.(vtotal)(Ia) = U_total.data(Ib);
        
        % Tidal Currents
        tab.(vdir_tidal)(Ia) = D_tide.data(Ib);
        tab.(vtidal)(Ia) = U_tide.data(Ib);
        
        % Residual Currents
        tab.(vdir_res)(Ia) = D_res.data(Ib);
        tab.(vres)(Ia) = U_res.data(Ib);
        
    end
    

    %% WRITE FILE
    
    % output file name
    output_fname = sprintf('ASOW_P1-DHI-MET-2022 Metocean data_Loc%0.3d_x.0.csv',loc); % per LOCATION
    
    % open file to write
    fid = fopen([output_dir,output_fname],'w');
    
    % add header
    head_info = {'<LOC>', sprintf('%0.3d',loc);
        '<LON>',sprintf('%0.3f',WL.xyz(1));
        '<LAT>',sprintf('%0.3f',WL.xyz(2));
        '<DEP>',sprintf('%0.1f m',WL.xyz(3));
        '<TODAY>',datestr(now(),'dd mmm yyyy')};
    
    % create header
    fprintf(fid,'DHI Water & Environment Inc.\n');
    fprintf(fid,'Project: ATLANTIC SHORES OFFSHORE WIND\n');
    fprintf(fid,'Date: %s\n', datestr(now(),'dd mmm yyyy'));
    
    fprintf(fid,'\nMetocean Time Series at Location %0.3d\n\n',loc);
    
    fprintf(fid,'Longitude: %0.3f\n',WL.xyz(1));
    fprintf(fid,'Latitude:  %0.3f\n',WL.xyz(2));
    fprintf(fid,'Depth:     %0.1f\n',WL.xyz(3));
    fprintf(fid,'\n');
    
    % add description
    for n = 1:height(params)
        fprintf(fid,'%-20s,%s\n',params.Parameter{n},params.Description{n});
    end
    
    fprintf(fid,'\n');
    fprintf(fid,'%25s', tab.Properties.VariableNames{1});
    for n = 2:length(tab.Properties.VariableNames)
        fprintf(fid,',%20s', tab.Properties.VariableNames{n});
    end
    fprintf(fid,'\n');
    for k = 1:height(tab)
        disp(k)
        fprintf(fid,'%25s',datestr(tab.Time(k),'yyyy/mm/dd HH:MM:SS'));
        dat = table2array(tab(k,2:end));
        dat(isnan(dat)) = -9999;
        fprintf(fid,',%20.2f', dat);
        fprintf(fid,'\n');
    end
    
    fclose(fid);
    
end

fclose all;