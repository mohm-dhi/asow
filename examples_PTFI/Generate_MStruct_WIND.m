%"""
%###############################################################################,
%# Created Date:    '2022-11-20'                                                #,
%# Author(s):       alka                           #,
%# Encoding:        UTF-8                                                       #,
%# Language:        MATLAB                                                      #,
%# ---------------------------------------------------------------------------- #,
%# Script for processing multiple stations with WL data from an HD model        #,
%#                                                                              #,
%#                                                                              #,
%# ---------------------------------------------------------------------------- #,
%# Copyright (c) (2022) DHI Water & Environment, Inc.                           #,
%################################################################################
%"""
%--------------------------------------------------------------------------
%01.- Preamble
%--------------------------------------------------------------------------
clc; clear all; close all; fclose all;

% Define main directory (either network or local)
MainDir = '//usden1-nas3/41806379/';
%MainDir = 'c:/Users/alka/Projects/41806379 - PTFI LNG Metocean/';

%% LOAD DATA
cfsr = readtable([MainDir,'02_Data/02_Processed/01_winds/CFSR_CFSv2/CFSR_136.56_-5.15_Wind.csv']);
cfsv2 = readtable([MainDir,'02_Data/02_Processed/01_winds/CFSR_CFSv2/CFSv2_136.64_-5.01_Wind.csv']);

time = datenum([cfsr.time; cfsv2.time]);
u = [cfsr.u_component_of_wind_height_above_ground; cfsv2.u_component_of_wind_height_above_ground];
v = [cfsr.v_component_of_wind_height_above_ground; cfsv2.v_component_of_wind_height_above_ground];

[ws,wd] = uv2spddir(u,v,'invert');

%--------------------------------------------------------------------------
%01- Winds
%--------------------------------------------------------------------------

data_fold = [MainDir,'02_Data/02_Processed/01_winds/CFSR_CFSv2/'];
dfs0_WS_name = ['CFSR_136.56_-5.15_Wind.csv'];

ttt_HDne =[datenum(1980,1,1) datenum(2021,1,1) 60 60];
year0 = 1980;
year1 = 2020;

xyz      =[136.6 -5.1 nan];

% General Definitions for the matlab structure
name            = 'Offshore';
bins_WS         = 0:2:20;     %inspection of data shows maximum values ~18m/s
bins_WD         = 0:22.5:360;

%========================Create structure =========================
legend_WS    = ['Wind Speed'];          %define legend for plots
legend_WD    = ['Wind Direction'];          %define legend for plots

WS       = m_structure(name,xyz,ttt_HDne,legend_WS, ...
    [time, ws]      ,'WS'      , 1 ,bins_WS);       %create m structure with all previously defined data

WD       = m_structure(name,xyz,ttt_HDne,legend_WS, ...
    [time, wd]      ,'WD'      , 1 ,bins_WS);       %create m structure with all previously defined data

WS.bins = bins_WS;
WD.bins = bins_WD;

WS.T  =[1,1.5,2,5,10,15,20,25,50,100];
WS.Tc =[1,1.5,2,5,10,15,20,25,50,100];

WS.xyz_str      = sprintf('(%0.1fE,%0.1fS)',abs(xyz(1)),abs(xyz(2)));
WS.ttt_str      = sprintf('(%d-%d)',year0,year1);
WS.ttt_str_long = sprintf('(%d-%d)',year0,year1);

WD.xyz_str      = sprintf('(%0.2fE,%0.2fS)',abs(xyz(1)),abs(xyz(2)));
WD.ttt_str      = sprintf('(%d-%d)',year0,year1);
WD.ttt_str_long = sprintf('(%d-%d)',year0,year1);

save([MainDir,'06_Results/M_STRUCTURES/' name '_WS.mat'], 'WS');           %Save Stations
save([MainDir,'06_Results/M_STRUCTURES/' name '_WD.mat'], 'WD');           %Save Stations
%m_timeseries_PTFI(WS,WD);



