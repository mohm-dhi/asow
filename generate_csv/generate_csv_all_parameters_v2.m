% generate deliverable for Atlantic Shores
% Wave spectral data: ASOW_P1-DHI-MET-yyyy Metocean data_LocXXXX_x.0.csv
% Where CONTRACTOR is the Contractor acronym, yyyy the year of delivery
% XXXX is the location name and
% x.0 is the revision status.
%
% ALKA, Dec 2022

clear all; close all; clc

FormatTable4Client = 1;

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
    
    % read header from mood
    header = readlines([sw_list.folder,'/',sw_list.name]);
    header = header(1:14);
    
    % ensure time matches and add to final table
    [~,Ia,Ib] = intersect(time,wave.datetime_ISO8601__UTC_);
    
    % translate MOOD header to final table
    for n = 1:height(mood2tab)
        tab.(mood2tab.Final{n})(Ia) = wave.(mood2tab.MOOD{n})(Ib);
    end
    
    %% WIND - 10m and hub height
    % load files
    load([wind_path,'ASOW',num2str(loc),'_WindSpd_10m.mat']);
    load([wind_path,'ASOW',num2str(loc),'_WindDir_10m.mat']);
    load([wind_path,'ASOW',num2str(loc),'_WindSpd_Hub.mat']);
    
    % ensure time matches and add to final table
    [~,Ia,Ib] = intersect(time,WindSpd_10m.datetime);
    
    tab.WindSpd_10m(Ia) = WindSpd_10m.data(Ib);
    tab.WindDir_10m(Ia) = WindDir_10m.data(Ib);
    tab.WindSpd_Hub(Ia) = WindSpd_Hub.data(Ib);
    tab.WindDir_Hub(Ia) = WindDir_10m.data(Ib);
    
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
    
    %% SAVE TABLE (easier QA)
    tab_fname = sprintf('BASE_TABLE_ASOW_P1-DHI-MET-2022 Metocean data_Loc%0.3d_x.0.csv',loc);
    writetable(tab,[output_dir,tab_fname])
    
    if FormatTable4Client
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
        fprintf(fid,'\nMetocean Time Series at Location %0.3d\n\n',loc);
        
        fprintf(fid,'Date Created: %s\n', datestr(now(),'dd mmm yyyy'));
        fprintf(fid,'DHI Project Manager: Danker Kolijn (dank@dhigroup.com)\n');
        fprintf(fid,'File Created by: Aline Kaji (alka@dhigroup.com)\n');
        fprintf(fid,'Point name: P%d\n',loc);
        fprintf(fid,'Longitude [deg.E]: %0.3f\n',WL.xyz(1));
        fprintf(fid,'Latitude[deg N]:  %0.3f\n',WL.xyz(2));
        fprintf(fid,'Projection: 4326\n');
        fprintf(fid,'Depth [mMSL]:     %0.1f\n',WL.xyz(3));
        fprintf(fid,'Start Date [UTC]: 1979-01-01 01:00:00\n');
        fprintf(fid,'End Date [UTC]→→: 2021-12-31 23:00:00\n');
        fprintf(fid,'Time Step [s]→→: 3600.0\n');
        
        fprintf(fid,'Dataset ID: US-EastCoast-SW, US-EastCoast-HD, CFSR-CFSv2\n')
        fprintf(fid,'This dataset is accompanied by a technical memorandum dated December 16, 2022 issued to ASOW (claire.dereux@atlanticshoreswind.com)\n')
        fprintf(fid,'The accompanying technical memorandum describes the various data sources and methods to derive the values contained in this file\n')
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
end

fclose all;