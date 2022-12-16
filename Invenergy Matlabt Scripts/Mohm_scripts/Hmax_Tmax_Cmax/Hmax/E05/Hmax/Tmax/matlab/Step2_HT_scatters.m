% Clean
clear all; close all; fclose all; clc

%% 0) User input
    ZC_file              = '\\DKCPH1-STOR2.DHI.DK\Projects\41806038_WS\RUGA\Hmax\E05\Hmax\Tmax\E_05_Zerocrossing_ED2f_deg_SW_{US_East_Coast,CFSR}_(1979-01-01_-_2020-12-31).mat';
    %SW_file              = '\\usden1-stor\Projects\41806038\Waves E05\US_East_Coast_Wave_Parameters_Integrated_MIKE_21_Spectral_Wave_Model_DHI_-72.71669_39.96928_NoHurricans.dfs0'; % For directional analysis
    SW_file              = '\\usden1-stor\Projects\41806038\Waves E05\US_East_Coast_Wave_Parameters_Integrated_MIKE_21_Spectral_Wave_Model_DHI_-72.71669_39.96928.dfs0'; % For directional analysis
    out_dir              = '\\DKCPH1-STOR2.DHI.DK\Projects\41806038_WS\RUGA\Hmax\E05\Hmax\Tmax';
    
    ttt                  = [datenum(1979,01,01,01,00,00) datenum(2020,12,31,23,00,00) 60 60];
    xyz                  = [-72.71669 39.96928 57];
    name                 = 'E 05';
    lgnd                 = 'SW_{US East Coast,CFSR}';
    
    bins_H               = [0:0.5:15];
    bins_T               = [0:0.5:27];
    bins_D               = [0:30:360];

%% 1) M structures for H, T and D
    H_SW                 = m_structure(name,xyz,ttt,lgnd,ZC_file,'H'        , 1,bins_H);
    T_SW                 = m_structure(name,xyz,ttt,lgnd,ZC_file,'T'        , 2,bins_T);
    
    [~,tmp]              = dfs02mat_dotnet(SW_file);
    
    % Remove empty registers
    Time                 = tmp.Time;
    MWD                  = tmp.Values(:,6);
%     ii                   = find(isnan(tmp.Values(:,6)));  
%     Time(ii)             = [];
%     MWD(ii)              = [];
   
    D_SW                 = m_structure(name,xyz,ttt,lgnd,[Time MWD ],'MWD_Total', 1,bins_D);
    
    H_SW.ascii           = '.xls';       
    T_SW.ascii           = '.xls';
    D_SW.ascii           = '.xls';
    
%% 2) Get scatter H-T
    root_dir             = pwd;
    cd(out_dir)
   % m_HT(H_SW,T_SW);
    m_HT(m_subseries(H_SW,'directional',D_SW),m_subseries(T_SW,'directional',D_SW));
    cd(root_dir)     
