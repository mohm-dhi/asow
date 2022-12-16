% Clean
clear all; close all; fclose all; clc

%% 0) User input
    SP_file              = '\\usden1-stor\Projects\41806038\Waves E05\SPEC\data\E05_SPEC.mat';
    out_dir              = '\\DKCPH1-STOR2.DHI.DK\Projects\41806038_WS\RUGA\Hmax\E05\Hmax\Tmax';
    
    ttt                  = [datenum(1979,01,01,01,00,00) datenum(2020,12,31,23,00,00) 60 60];
    xyz                  = [-72.71669 39.96928 57];
    name                 = 'E 05';
    lgnd                 = 'SW_{US East Coast,CFSR}';
    
    bins_ED2             = [0:10 :1000];   
    
%% 1) Step 01: Get the zerocrossing file out of spectra
    ED2_SW               = m_structure(name,xyz,ttt,lgnd,SP_file,'ED2f_deg',1,bins_ED2); 
    load(SP_file,'-mat'); 
    ED2_SW.data          = X.data; 
    ED2_SW.time          = X.time; 
    ED2_SW.dfs_unit      = 'm^2*s/deg';
    ED2_SW.xaxisUnit     = X.xaxisUnit; 
    ED2_SW.yaxisUnit     = X.yaxisUnit; 
    ED2_SW.dataStructure = X.dataStructure; 
    ED2_SW.dfs_type      = X.dfs_type;
    ED2_SW.dfs_name      = X.dfs_name;     
    
    % Run m_zerocrossing
    root_dir             = pwd;
    cd(out_dir)
    m_zerocrossing(ED2_SW);
    cd(root_dir)      