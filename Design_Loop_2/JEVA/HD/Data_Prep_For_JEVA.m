%"""
%################################################################################,
%# Created Date:    '2023-04-15'                                                #,
%# Author(s):       Seyed Abbas Jazaeri - abja@dhigroup.com                     #,
%# Encoding:        UTF-8                                                       #,
%# Language:        MATLAB R2022b                                               #,
%# ---------------------------------------------------------------------------- #,
%# This script prepares the input files required for JEVA                       #,
%# IMPORTANT NOTICE: Please note that to have consistant MATLAB files with      #,
%# the input file required for the JEVA file, some minor modifications have     #,
%# implemented in "m_tide" function. For more information, take a look at       #,
%# the README file or contact abja@dhigroup.com.                                #,
%# ---------------------------------------------------------------------------- #,
%# Copyright (c) (2022) DHI Water & Environment, Inc.                           #,
%################################################################################
%"""

close all
clear
clc

%=============================================================================================================
% Definition of the required input data 
%=============================================================================================================

% fname = '0D_model_cable_corr_cmp_1980_2022.dfs0';
% stations = {'CC-01', 'CC-02', 'CC-03', 'CC-04', 'CC-05', 'CC-06'};
% lons = [-73.995, -73.992, -74.030, -74.027, -74.441, -74.403];
% lats = [40.224, 40.224, 40.116, 40.113, 39.345, 39.331];
% depths = [-5.2, -10.7, -5.2, -10.7, -5.2, -10.8];
% ttt = [datenum([1979 12 15 00 00 00]) datenum([2023 01 01 00 00 00]) 60];

fname = '0D_model_ocs_a0499_a0545_cmp_1980_2022.dfs0';
stations = {'LA-01', 'LA-02', 'LA-03', 'LA-04', 'LA-05', 'LA-06', 'LA-07', 'LA-08', 'LA-09', 'LA-10', 'LA-11', 'LA-12', 'OSSA' , 'OSSB'};
lons =     [-73.960, -74.034, -73.953, -73.946, -74.100, -74.021, -73.953, -74.236, -74.052, -74.087, -74.026, -74.069, -74.166, -74.136];
lats =     [39.619 , 39.559 , 39.484 , 39.367 , 39.364 , 39.340 , 39.265 , 39.279 , 39.201 , 39.147 , 39.424 , 39.283 , 39.288 , 39.224];
depths =   [-24.8  , -23.6  , -25.5  , -27.6  , -22.8  , -27.0  , -31.5  , -19.2  , -22.7  , -36.5  , -25.4  , -30.3  , -26.3  , -23.4];
ttt = [datenum([1979 12 15 00 00 00]) datenum([2023 01 01 00 00 00]) 60];

% fname = '0D_model_ocs_ocs_a0541_cmp_1980_2022.dfs0';
% stations = {'LA-21', 'LA-22', 'LA-23', 'LA-24', 'LA-25'};
% lons = [-73.453, -73.510, -73.650, -73.618, -73.579];
% lats = [39.439, 39.334, 39.193, 39.485, 39.362];
% depths = [-34.2, -50.1, -43.2, -35.5, -46.4];
% ttt = [datenum([1979 12 15 00 00 00]) datenum([2023 01 01 00 00 00]) 60];

% The definition of the bins for each variable
binsWL = -4:0.5:4;
binsUV = -1:0.5:1;
binsCS = 0:0.1:1.5;
binsCD = 0:30:330;

currd = pwd;
mkdir(fullfile(currd, "m_tide_outputs"))
mkdir(fullfile(currd, "HD"))
cd(fullfile(currd, "m_tide_outputs"))

for i = 1:length(stations)
    mkdir(fullfile(currd, "m_tide_outputs", stations{i}));
    cd(fullfile(currd, "m_tide_outputs", stations{i}));

%=============================================================================================================
% Creating the m_structure variables based on the input data
%=============================================================================================================
    xyz = [lons(i) lats(i) depths(i)];

    WL = m_structure(stations{i} , xyz , ttt , 'HD_{US-EC}' , fullfile(currd, fname) ,'WLtot' , i, binsWL);
    US = m_structure(stations{i} , xyz , ttt , 'HD_{US-EC}' , fullfile(currd, fname) ,'U' , 1*length(stations)+i, binsUV);
    VS = m_structure(stations{i} , xyz , ttt , 'HD_{US-EC}' , fullfile(currd, fname) ,'V' , 2*length(stations)+i, binsUV);

%=============================================================================================================
% De-tiding the water level and current data
%=============================================================================================================

    WL_Components = m_tide(WL);
    WLtot = WL_Components.WL;
    WLtide = WL_Components.WLTde;
    WLres = WL_Components.WLRes;
        
    disp(strcat(stations{i}, " water level de-tided!"))

    Current_Components = m_tide(US, VS);

    disp(strcat(stations{i}, " current speed and direction de-tided!"))

    CStot = Current_Components.CS;
    CStide = Current_Components.CSTde;
    CSres = Current_Components.CSRes;
    CDtot = Current_Components.CD;
    CDtide = Current_Components.CDTde;
    CDres = Current_Components.CDRes;

%=============================================================================================================
% Prepare mat file for JEVA
%=============================================================================================================
%     The data are selected with every other one pattern because the JEVA
%     code takes hourly data but the provided dfs0 files are half hourly
    data = [WLtot.data(1:2:end) WLtide.data(1:2:end) WLres.data(1:2:end) CStot.data(1:2:end) CStide.data(1:2:end) CSres.data(1:2:end) CDtot.data(1:2:end) CDtide.data(1:2:end) CDres.data(1:2:end)];
    item = {WLtot.item WLtide.item WLres.item CStot.item CStide.item CSres.item CDtot.item CDtide.item CDres.item};

    metadata.Project_Name         = 'ASOW - Design Loop 2';     % Check this value!
    metadata.Produced_by          = 'DHI A/S';                  % Check this value!
    metadata.Accompanying_report  = 'Report.docx';              % Check this value!
    metadata.Data_prepared_by     = 'ABJA';
    metadata.Data_checked_by      = 'ALKA';                     % Check this value!
    metadata.Date_of_preparation  =  datestr(today("datetime"), "yyyy-mm-dd");
    metadata.Model_Version        = 'HD_{US-EC}';               % Check this value!
    metadata.Point_Name           = stations{i};
    metadata.Projection           = 'Geographical Long/Lat';
    metadata.Longitude            = lons(i);
    metadata.Latitude             = lats(i);
    metadata.Modelled_Water_Depth = depths(i);
    metadata.Time_Zone            = 'UTC';
    metadata.Mean_Mesh_Resolution = '400 m';                    % Check this value!

    time = WLtot.time(1:2:end);
%     The units of the directions are hard-coded to be consistant with the
%     JEVA sample .mat file
    unit = {WLtot.unit WLtide.unit WLres.unit CStot.unit CStide.unit CSres.unit 'Deg. N. (going-to)' 'Deg. N. (going-to)' 'Deg. N. (going-to)'};

    cd(fullfile(currd, "HD"));
    save(strcat('HD_', stations{i}, '.mat'), "data", "item", "metadata", "time", "unit");

    disp(strcat(stations{i}, " JEVA MATLAB file created!"))
%=============================================================================================================

end
