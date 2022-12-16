%"""
%###############################################################################,
%# Created Date:    '2022-10-22'                                                #,
%# Author(s):       Pablo Cortes - pcsa@dhigroup.com                            #,
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
% addpath('c:\Users\alka\Documents\MATLAB\DHI-MATLAB-Toolbox-master\mbin\')
addpath('DHI-MATLAB-Toolbox\mbin\')
addpath('m_tools_local\')
addpath('potlab_v2\')

NET.addAssembly('DHI.Mike.Install');
import DHI.Mike.Install.*;
DHI.Mike.Install.MikeImport.SetupLatest({DHI.Mike.Install.MikeProducts.MikeCore});

Stations=readtable('Stations_Mstructs_2022_11_22.xlsx');
currd = pwd;   
outdir = '../06_Results/M_STRUCTURES/';
load('hycom_signal/HYCOM_HF_SIGNAL.mat');

%--------------------------------------------------------------------------
%HD Matlab structure creation
%--------------------------------------------------------------------------    

year0 = 2000;
year1 = 2020;

proc_list=[2,4,24,25]; %N molf, S molf, Fluor S Molf, 2MB, FSRU, 

    bins_WL = -2.75:0.25:2.75;   
    bins_CS = 0:0.2:1.4; 
    bins_CD = 0:22.5:360; 

    data_fold_wl = '../../../06_Results/02_HD/41806379_Hindcast_extractions/WL/';
    data_fold_csd = '../../../06_Results/02_HD/41806379_Hindcast_extractions/CSD/';

    dfs0_HD_name = 'mfm_ext_conc_2022_11_22.dfs0';
%     dfs0_SW_name = 'mfm_ext_conc_2022_11_22.dfs0';

    for i=1:size(proc_list,2)     
        ttt_HDne =[ datenum(Stations.StartDate(proc_list(i))) datenum(Stations.EndDate(proc_list(i))) ...
                    Stations.Delta__min_(proc_list(i))*ones(1,3)];

        xyz      =[Stations.Longitude(proc_list(i)) Stations.Latitude(proc_list(i)) Stations.ModelWaterDepth(proc_list(i))];

        name     = [Stations.Label{proc_list(i)}]; 
        %========================Create structure, WL =====================
        legend_WL    = ['Water Level'];                                                             %define legend for plots
        WL       = m_structure(name,xyz,ttt_HDne,legend_WL, ...      
            [data_fold_wl dfs0_HD_name ]      ,'WL'      , proc_list(i)-1 ,bins_WL);                   %-1 since station 1 is wind offshore
        
%        WL.bins=WL.bins';                                                                     %fix bin from column to row (incompatibility with some fnctions) 
        WL.ttt_str     =''; WL.ttt_str_long='';                                               %remove time string for plots            
        WL.data= WL.data-HD_HYCOM_FLT.data;                                             %subtract hycom filtered signal to avoid double counting of surge
        
        WL.ttt_str      = sprintf('(%d-%d)',year0,year1);

        
        WL.T  =[1,1.5,2,5,10,15,20,25,50,100]; 
        WL.Tc =[1,1.5,2,5,10,15,20,25,50,100]; 
        
        save([outdir Stations.Label{proc_list(i)} '_WL.mat'], 'WL');            %Save Stations

        %---------------------Create tidal structure-----------------------
        X=m_resample(WL,'dt_int',60,'dt_avg',10,'maxgap',10);            %lighter structure for faster tidal calculations
        [Utide_tide,Utide_resi,WL_Tidal_Coef_model]=m_tide_utide_local(X);                   %Residual and tide every 1h because X is in 1h       


        
        WL_Tidal_Pred_model=m_resample(Utide_tide,'dt_int',10,'dt_avg',60,'maxgap',60); %Resample tide every 10 min        
        WL_Tidal_Pred_model.name=[Stations.Label{proc_list(i)} ' Tidal Prediction'];
        
        save([outdir Stations.Label{proc_list(i)} '_WL_Tidal_Pred_model.mat'], 'WL_Tidal_Pred_model'); %Save utide tidal prediction          
        save([outdir Stations.Label{proc_list(i)} '_WL_Tidal_Coef_model.mat'], 'WL_Tidal_Coef_model'); %Save utide tidal prediction                      
        
         %------------------Create Skew Residual Structures----------------          
        [WL_Res_Skew_HT,WL_Res_Skew_LT]=m_skew_residual(WL,WL_Tidal_Pred_model);        
        WL_Res_Skew_HT.T=WL.T; WL_Res_Skew_HT.Tc=WL.Tc;
        WL_Res_Skew_LT.T=WL.T; WL_Res_Skew_LT.Tc=WL.Tc;
        WL_Res_Skew_HT_int.bins=-1:0.1:1;
        WL_Res_Skew_LT_int.bins=-1:0.1:1;
        save([outdir Stations.Label{proc_list(i)} '_WL_Res_Skew_HT.mat'], 'WL_Res_Skew_HT'); %Save utide tidal prediction          
        save([outdir Stations.Label{proc_list(i)} '_WL_Res_Skew_LT.mat'], 'WL_Res_Skew_LT'); %Save utide tidal prediction                  
        %========================Create structure, CSD ====================
        legend_CS    = ['Current Speed'];                                                           %define legend for plots
        CS     = m_structure(name,xyz,ttt_HDne,legend_CS, ...
                     [data_fold_csd dfs0_HD_name ]      ,'CS'      ,(proc_list(i)-1) ,bins_CS);      %create m structure with all previously defined data          

        legend_CD    = ['Current Direction, From'];                                                 %define legend for plots
        CD     = m_structure(name,xyz,ttt_HDne,legend_CD, ...
                     [data_fold_csd dfs0_HD_name ]      ,'CD'      ,24+(proc_list(i)-1)  ,bins_CD);     %create m structure with all previously defined data   

        CS.T  =[1,1.5,2,5,10,15,20,25,50,100]; 
        CS.Tc =[1,1.5,2,5,10,15,20,25,50,100];                  
                 
        CS.ttt_str     =sprintf('(%d-%d)',year0,year1);
        CS.ttt_str_long=sprintf('(%d-%d)',year0,year1);                            %remove time string for plots   
        CD.ttt_str     =sprintf('(%d-%d)',year0,year1);
        CD.ttt_str_long=sprintf('(%d-%d)',year0,year1);                            %remove time string for plots   
        save([outdir Stations.Label{proc_list(i)} '_CS.mat'], 'CS');       %Save Stations
        save([outdir Stations.Label{proc_list(i)} '_CD.mat'], 'CD');       %Save Stations
        
        %------------------------------------------------------------------
        %Estimate surface and bottom current speeds.
        %------------------------------------------------------------------
        h=abs(Stations.ModelWaterDepth(proc_list(i))); %local depth. Check if use bathy or model bathy

        CS_Surface = CS; CS_Bottom  = CS;
        
        %Soulsby 1997.
        alph= 7; beta= 0.32; %Shelf Sea.
        z   = 1;             %use fixed 1m height above bed.
        
        CS_Surface.data = 1.07                 *CS.data; % Z above h/2, h equals local depth
        CS_Bottom.data  = (z/(beta*h))^(1/alph)*CS.data;        
        
        CS_Surface.legend= 'Surface Current Speed';
        CS_Bottom.legend=  'Bottom  Current Speed';
        
        save([outdir Stations.Label{proc_list(i)} '_CS_Surface.mat'], 'CS_Surface');       %Save Stations
        save([outdir Stations.Label{proc_list(i)} '_CS_Bottom.mat'], 'CS_Bottom');       %Save Stations        
        
        %========================TIMESERIES (OPTIONAL) ==================

        m_timeseries_PTFI(WL);
        m_timeseries_PTFI(WL_Tidal_Pred_model)
        m_timeseries_PTFI(CS,CD);                       
        m_timeseries_PTFI(CS,CS_Surface,CS_Bottom); 
        cd(currd)
    end