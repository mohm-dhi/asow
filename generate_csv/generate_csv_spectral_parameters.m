% generate deliverable for Atlantic Shores
% Wave spectral data: ASOW_P1-DHI-MET-yyyy Metocean data_LocXXXX_x.0.csv
% Where CONTRACTOR is the Contractor acronym, yyyy the year of delivery
% XXXX is the location name and
% x.0 is the revision status.
%
% ALKA, Dec 2022

clear all; close all; clc

%% SETTINGS
% time series vector (to be uniform for all parameters)
time = (datetime(1979,1,15):hours(1):datetime(2021,12,31,23,0,0))';

% output directory
output_dir = '\\USDEN1-STOR.DHI.DK\Projects\41806529\07_Timeseries_CSV_Deliverables\';

% template
template_fname = '\\USDEN1-STOR.DHI.DK\Projects\41806529\07_Timeseries_CSV_Deliverables\TEMPLATE_TimeSeries';

% list of parameters
params_fname = '\\USDEN1-STOR.DHI.DK\Projects\41806529\07_Timeseries_CSV_Deliverables\ListParameters.xlsx';

% wave parameters from 2D spectra
spec2d_Dir = '\\USDEN1-STOR.DHI.DK\Projects\41806529\03_PROCESSED_MOOD_data\02-OutputLocations\SpectralAnalysis\';

%% READ TEMPLATE
template = readlines(template_fname);

% read overall parameters table
params = readtable(params_fname,'Sheet','SpectralParams');

%% LOOP THROUGH OUTPUT LOCATIONS
for loc = 8
    
    %% LOAD 2D SPECTRA PARAMETERS
    load([spec2d_Dir,sprintf('Spectra_Integral_Parameters_Loc%0.3d.mat',loc)],'S');
    
    %% INITIALIZE OUTPUT TABLE
    tab = array2table(nan(length(time),height(params)),...
        'VariableNames',params.Parameter);
    tab.Time = time;
    
    % ensure time matches and add to final table
    [~,Ia,Ib] = intersect(time,S.Time);
    
    % add data to table
    for n = 1:height(params)
        varname = params.Parameter{n};
        if isfield(S,varname)
            tab.(varname)(Ia) = S.(varname)(Ib);
        end
    end
    
    %% WRITE FILE
    
    % output file name
    output_fname = sprintf('ASOW_P1-DHI-MET-2022 Integral spectral data_Loc%0.3d_x.0.csv',loc); % per LOCATION
    
    % open file to write
    fid = fopen([output_dir,output_fname],'w');

    % create header
    fprintf(fid,'DHI Water & Environment Inc.\n');
    fprintf(fid,'Project: ATLANTIC SHORES OFFSHORE WIND\n');
    fprintf(fid,'Date: %s\n', datestr(now(),'dd mmm yyyy'));
    
    fprintf(fid,'\nMetocean Time Series at Location %0.3d\n\n',loc);
    
    fprintf(fid,'Longitude: %0.3f\n',S.xyz(1));
    fprintf(fid,'Latitude:  %0.3f\n',S.xyz(2));
    fprintf(fid,'Depth:     %0.1f\n',S.xyz(3));
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