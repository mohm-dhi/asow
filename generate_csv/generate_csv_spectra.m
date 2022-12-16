% generate deliverable for Atlantic Shores
% Wave spectral data: ASOW_P1-CONTRACTOR-MET-yyyy Wave spectral data_LocXXXX_x.0.csv
% Where CONTRACTOR is the Contractor acronym, yyyy the year of delivery
% XXXX is the location name and
% x.0 is the revision status.
%
% ALKA, Dec 2022

clear all; close all; clc

% load data
load('\\USDEN1-STOR.DHI.DK\Projects\41806529\02_RAW_MOOD_data\2D Spectra data\Spec0.1_1979-2021_-73.9W39.3N.mat');

% output file name
fname = 'ASOW_P1-DHI-MET-2022 Wave spectral data_Loc000_%0.4d_x.0.csv'; % per year

% output directory
output_dir = '\\USDEN1-STOR.DHI.DK\Projects\41806529\07_Timeseries_CSV_Deliverables\';

% depth at location
% info from ELME
Depth = '33.4 mMSL';

for year = 1979:2021
    
    disp(year)
    
    % select only data for this year
    Xyear = X;
    mask = (X.time >= datenum(year,1,1)) & (X.time < datenum(year+1,1,1));
    Xyear.time = X.time(mask);
    Xyear.data = X.data(mask,:,:);
    
    % add year to filename
    fname_year = sprintf(fname,year);
    
    %% START OUTPUT FILE
    fid = fopen([output_dir,fname_year],'w');
    
    %% WRITE GENERAL HEADER
    fprintf(fid,'%% DHI Water & Environment Inc.\n');
    fprintf(fid,'%% Project: %s\n', Xyear.Project_Name);
    fprintf(fid,'%% Date: %s\n', datestr(now(),'dd mmm yyyy'));
    fprintf(fid,'%% Parameter: 2D Wave Spectra\n');
    fprintf(fid,'%% Model: U.S. East Coast SW %s\n', Xyear.Model_Version);
    fprintf(fid,'%% Longitude: %s\n',Xyear.Easting);
    fprintf(fid,'%% Latitude: %s\n',Xyear.Northing);
    fprintf(fid,'%% Depth: %s\n\n',Depth);
    
    fprintf(fid,'FREQUENCIES [Hz] - rows\n');
    fprintf(fid,'%0.5f\n',Xyear.xaxis);
    fprintf(fid,'DIRECTIONS [degN] - columns\n');
    fprintf(fid,'%0.5f\n',Xyear.yaxis);
    
    %% WRITE
    for t = 1:length(Xyear.time)
        
        if ~mod(t,400)
            disp(datestr(Xyear.time(t)))
        end
        
        fprintf(fid,'\nTIME [UTC]: %s\n', datestr(Xyear.time(t),'yyyy-mm-dd HH:MM'));
        fprintf(fid,'ENERGY [m2/Hz/deg]\n'); % TO ADD - CONFIRM UNITS
        for k = 1:size(Xyear.data,2)
            for m = 1:size(Xyear.data,3)
                value = squeeze(Xyear.data(t,k,m));
                if m > 1
                    fprintf(fid,',');
                end
                fprintf(fid,'%0.6f',value);
            end
            fprintf(fid,'\n');
        end
    end
    
    fclose(fid);
    
end