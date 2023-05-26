close all;
clear;
clc;

%NET.addAssembly('DHI.Mike.Install');
%import DHI.Mike.Install.*;
%DHI.Mike.Install.MikeImport.SetupLatest({DHI.Mike.Install.MikeProducts.MikeCore});

config_file_name = 'config.csv';

opts = detectImportOptions(config_file_name);
opts.Delimiter = {','};
opts.VariableTypes{2} = 'char'; % make the stations column data type to char

config = readtable(config_file_name, opts);
config = config(string(config.disable)=="", :);

stations = config.stations;
types = config.type;
locs = [config.lon config.lat config.depth];
obs_path = config.obs_path;
model_path = config.model_path;
% obs_column = str2num(char(string(config.obs_column)));
% model_column = str2num(char(string(config.model_column)));
currd = pwd;


for i = 1:length(stations)
    directory = [stations{i},'_' ,types{i}];
    mkdir(fullfile(currd,directory));
    cd(fullfile(currd,directory));
    % General
    name            = stations{i};
    type            = types{i};
    xyz             = locs(i,:);
    ttt_obs        = [datenum([2021 05 01 00 00 00]) datenum([2022 04 01 00 00 00]) 10 10];
    ttt_model    	= [datenum([2021 05 01 00 00 00]) datenum([2022 04 01 00 00 00]) 60 120];
    ttt_model_hd    	= [datenum([2021 05 01 00 00 00]) datenum([2022 04 01 00 00 00]) 30 30];
    ttt_model_sw    	= [datenum([2021 05 01 00 00 00]) datenum([2022 04 01 00 00 00]) 60 180];
    if contains(type, 'wl')
%        ttt_obs     = [datenum([2002 3 1 00 00 00]) datenum([2022 12 01 00 00 00]) 6];
%        ttt_model    	= [datenum([2002 3 1 00 00 00]) datenum([2022 12 01 00 00 00]) 30];
       
       ttt_obs     = [datenum([2021 9 1 00 00 00]) datenum([2022 1 15 00 00 00]) 10 10];
       ttt_model    	= [datenum([2021 9 1 00 00 00]) datenum([2022 1 15 00 00 00]) 30 10];
    end
    % Bins
    bins_H         = 0:6;
    bins_T         = 0:2:24;
    bins_D         = 0:30:330;
    bins_WS        = 0:3:33;
    bins_CS        = 0:0.1:0.9;
    bins_WL        = -1.6:0.2:1.6;
    % U V
    columns_obs	    = str2num(char(string(config.obs_column(i))));
    columns_model	= str2num(char(string(config.model_column(i))));
    % Legends
    legend_obs_wind     = 'ASOW FLiDAR';
    legend_obs_hd     = 'Observations';
    legend_obs_sw     = 'Observations';
    legend_model_wind     = 'CFSR';
    legend_model_hd     = 'HD_{US-EC}';
    legend_model_sw     = 'SW_{US-EC}';



    % Files - Processed
    fname_model	=	model_path{i,1};
    fname_obs	=	obs_path{i,1};

    save cmp_data.mat name xyz type ttt_obs ttt_model ttt_model_hd ttt_model_sw bins_H bins_D columns_obs bins_WL columns_model bins_WS bins_T bins_CS legend_obs_wind legend_obs_hd legend_obs_sw legend_model_wind legend_model_hd legend_model_sw fname_model fname_obs
    m_file = ['t_' replace(stations{i}, ' ', '_') '_' types{i}];
    m_file = replace(m_file, '-', '_');
    copyfile(fullfile(currd,'Template.m'),[m_file,'.m']);
    eval(m_file)
    cd(currd)
end

