% Clean
close all; clear all; clc;

% Add paths
addpath(genpath('C:\Ruben_work\Matlab_Toolboxes'))

% HEADS UP! HARCODED LINES 57-61, 83-87 and 112-116

%% 0) User input
    modes_dir    = '\\DKCPH1-STOR2.DHI.DK\Projects\41806038_WS\RUGA\Hmax\E05\Modes';
    out_pth_Hmax = '\\DKCPH1-STOR2.DHI.DK\Projects\41806038_WS\RUGA\Hmax\E05\Hmax';
    out_pth_Cmax = '\\DKCPH1-STOR2.DHI.DK\Projects\41806038_WS\RUGA\Hmax\E05\Cmax';
    
    % EVA Hmax / Cmax
    EVAOpt.EVdist           = 3;                                           % 2. Trunc. Weibull 3. 2-p Weibull 4. Gumbel 5. Exponential
    EVAOpt.EVtype           = 'AAP';                                       % Annual Average Peak, quite similar to POT method
    EVAOpt.EVcrit           = 3;                                           % Number of events per year
    EVAOpt.N_BootStrap      = 500;                                         % Bootstraping to produce 95% confidence limit
    EVAOpt.estimationmethod = 'ML';                                        % Least square method to fit extremes onto the distribution
    EVAOpt.plotflag         = 2;
    EVAOpt.EVadju           = 0;
    EVAOpt.T                = [1 5 10 50 100 500 1000 10000];	           % Return period
    EVAOpt.intereventtime   = 36;	                                       % Inter event time between two extreme events (36 hours)
    EVAOpt.intereventlevel  = 0.7;
    EVAOpt.ndec             = 1;
    EVAOpt.ConfLimits       = [0.025 0.975];                               % 95% confidence
        
%% 1) Get files in modes_dir
    [Nfiles_H,file_list_H]     = ListFilesInDir_v2(modes_dir,'_H_Directional.mat');
    [Nfiles_MSL,file_list_MSL] = ListFilesInDir_v2(modes_dir,'_MSL_Directional.mat');
    [Nfiles_SWL,file_list_SWL] = ListFilesInDir_v2(modes_dir,'_SWL_Directional.mat');
    
    % Sort files  
    file_list_H                = natsortfiles(file_list_H);
    file_list_MSL              = natsortfiles(file_list_MSL);
    file_list_SWL              = natsortfiles(file_list_SWL);    

%% 2) Compute Cmax

    % 2.0) Root dir
         root_dir = pwd;

    % 2.1) Crests w.r.t. MSL 
    
         % Move to dir
         cd(out_pth_Cmax)
    
         % Start waitbar
         h = waitbar(0,'','Color',[0.5 0.6 0.7],'Name','Processing points for Cmax MSL');   
    
         for i = 1:Nfiles_MSL
             
             load(fullfile(modes_dir,file_list_MSL{i}))
             %EVAOpt.xtick = 6:0.5:11; 
             EVA_Cmax_MSL = m_extreme(modedata,EVAOpt);
             %outname      = strcat('Location_',num2str(i),'_Cmax_MSL.mat');
             if i ==1
                outname    = 'NY_Bight_Forristall_Cmax_MSL.mat';
             else
                outname    = 'NY_Bight_Glukhovskiy_Cmax_MSL.mat';
             end
             save(fullfile(out_pth_Cmax,outname),'EVA_Cmax_MSL')
             
             % Update waitbar
             waitbar((i/Nfiles_MSL),h,strcat('Progress: ',{' '},num2str((i/Nfiles_MSL)*100,3),'%',' ','(',num2str(i),' of ',num2str(Nfiles_MSL),')'))                  
           
         end
         
         % Close waitbar
         close(h)

    % 2.2) Crests w.r.t. SWL

         % Start waitbar
         h = waitbar(0,'','Color',[0.5 0.6 0.7],'Name','Processing points for Cmax SWL');       
    
         for i = 1:Nfiles_SWL
             
             load(fullfile(modes_dir,file_list_SWL{i}))
             %EVAOpt.xtick = 6:0.5:11; 
             EVA_Cmax_SWL = m_extreme(modedata,EVAOpt);
             %outname      = strcat('Location_',num2str(i),'_Cmax_SWL.mat');
             if i ==1
                outname    = 'NY_Bight_Forristall_Cmax_SWL';
             else
                outname    = 'NY_Bight_Glukhovskiy_Cmax_SWL';
             end             
             save(fullfile(out_pth_Cmax,outname),'EVA_Cmax_SWL')
             
             % Update waitbar
             waitbar((i/Nfiles_SWL),h,strcat('Progress: ',{' '},num2str((i/Nfiles_SWL)*100,3),'%',' ','(',num2str(i),' of ',num2str(Nfiles_SWL),')'))                 
           
         end
         
         % Close waitbar
         close(h)
         
    % 2.3) Hmax
        
         % Move to dir
         cd(out_pth_Hmax)        
    
         % Start waitbar
         h = waitbar(0,'','Color',[0.5 0.6 0.7],'Name','Processing points for Hmax');       
    
         for i = 1:Nfiles_H
             
             load(fullfile(modes_dir,file_list_H{i}))
             %EVAOpt.xtick = 6:0.5:11; 
             EVA_Hmax = m_extreme(modedata,EVAOpt);
             %outname  = strcat('Location_',num2str(i),'_Hmax.mat');
             if i ==1
                outname    = 'NY_Bight_Forristall_Hmax.mat';
             else
                outname    = 'NY_Bight_Glukhovskiy_Hmax.mat';
             end             
             save(fullfile(out_pth_Hmax,outname),'EVA_Hmax')
             
             % Update waitbar
             waitbar((i/Nfiles_H),h,strcat('Progress: ',{' '},num2str((i/Nfiles_H)*100,3),'%',' ','(',num2str(i),' of ',num2str(Nfiles_H),')'))                 
           
         end         

         % Close waitbar
         close(h)    
         
         % Move to Root dir
         cd(root_dir)    
    