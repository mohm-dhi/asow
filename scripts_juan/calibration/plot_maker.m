%"""
%###############################################################################,
%# Created Date:    '2022-06-21'                                                #,
%# Author(s):       Mohammad Madani - mohm@dhigroup.com                         #,
%#                  Benjamin Hernandez Alfaro - beha@dhigroup.com               #,
%# Encoding:        UTF-8                                                       #,
%# Language:        MATLAB 3.9                                                  #,
%# ---------------------------------------------------------------------------- #,
%# This script reads the configuration file and prepares the variables          #,
%#                                                                              #,
%#                                                                              #,
%# ---------------------------------------------------------------------------- #,
%# Copyright (c) (2022) DHI Water & Environment, Inc.                           #,
%################################################################################
%"""

clc; clear; close all;

addpath(genpath('C:\Users\jngz\MATLAB\Projects\potlab\src'));
addpath(genpath('C:\Users\jngz\MATLAB\Projects\potlab\res'));

NET.addAssembly('DHI.Mike.Install');
import DHI.Mike.Install.*;
DHI.Mike.Install.MikeImport.SetupLatest({DHI.Mike.Install.MikeProducts.MikeCore});

config_file_name = 'config.csv';

opts = detectImportOptions(config_file_name);
opts.VariableTypes{2} = 'char'; % make the stations column data type to char
config = readtable(config_file_name, opts);
config = config(string(config.disable)=="", :);

stations = config.stations;
locs = [config.lon config.lat config.depth];
obs_path = config.obs_path;
model_path = config.model_path;
obs_column = str2num(char(string(config.obs_column)));
model_column = str2num(char(string(config.model_column)));
currd = pwd;

for i = 1:length(stations)
    mkdir(fullfile(currd,stations{i}));
    cd(fullfile(currd,stations{i}));

    % General
    name            = stations{i};
    xyz             = locs(i,:);
    ttt_obs        = [datenum([2019 12 29 15 20 00]) datenum([2020 06 26 00 00 00]) 10];
    ttt_model    	= [datenum([1979 01 01 00 00 00]) datenum([2021 12 31 00 00 00]) 60];
    
    % Bins
    bins_H         = 0:1:22;
    bins_T         = 0:1:24;
    bins_D         = 0:10:360;
    bins_WS         = 0:3:33;
    bins_CS         = 0:0.1:1;
    bins_hm0        = 0:1:10;

    % U V
    colums_obs	    = obs_column(i,:);
    colums_model	= model_column(i,:);

    % Legends
    legend_obs     = 'Obs';
    legend_model     = 'Model';
        
    % Files - Processed
    fname_model	=	model_path{i,1};
    fname_obs	=	obs_path{i,1};
    
    save cmp_data.mat name xyz ttt_obs ttt_model bins_H bins_D colums_obs colums_model bins_CS bins_WS bins_T bins_hm0 legend_obs legend_model fname_model fname_obs
    
    copyfile(fullfile(currd,'template.m'),['t_',stations{i},'.m']);
	eval(['t_',stations{i}])
    
    cd(currd)
end
