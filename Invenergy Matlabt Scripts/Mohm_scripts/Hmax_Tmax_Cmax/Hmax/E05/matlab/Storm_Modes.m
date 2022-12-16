% Clean
close all; clear all; clc;

%% 0) User input
    file_SW     = '\\usden1-stor\Projects\41806038\Waves E05\US_East_Coast_Wave_Parameters_Integrated_MIKE_21_Spectral_Wave_Model_DHI_-72.71669_39.96928_NoHurricans.dfs0';
    file_HD     = '\\usden1-stor\Projects\41806038\HD E05\US_East_Coast_Water_Level_and_Current_2D_MIKE_21_Hydrodynamic_Model_HD_DHI_-72.71669_39.96928_NoHurricans.dfs0';
    out_pth     = '\\DKCPH1-STOR2.DHI.DK\Projects\41806038_WS\RUGA\Hmax\E05\Modes';

    Names       = {'E 05'};
    Lat         = [39.96928];
    Lon         = [-72.71669];
    Depth       = [57];      
    legend_SW   = 'SW_{US East Coast,CFSR}';
    legend_HD   = 'HD_{US East Coast,CFSR}';
    
    bins_Hm0    = [0:0.2:5]; 
    bins_Tp     = [0:1:30];
    bins_T01    = [0:1:20];    
    bins_T02    = [0:1:20];
    
    bins_WL     = [-2.0:0.25:2];  
    bins_PWD    = 0:30:360;
    bins_MWD    = 0:30:360;
    
%% 1) Load files
    [~,tmp_SW]  = dfs02mat_dotnet(file_SW);
    [~,tmp_HD]  = dfs02mat_dotnet(file_HD);   

%% 2) Get relevant info
    Time_SW    = tmp_SW.Time;
    Hm0        = tmp_SW.Values(:,1);
    Tp         = tmp_SW.Values(:,2);
    T01        = tmp_SW.Values(:,3);
    T02        = tmp_SW.Values(:,4);
    PWD        = tmp_SW.Values(:,5);
    MWD        = tmp_SW.Values(:,6);

    Time_HD    = tmp_HD.Time;
    WL         = tmp_HD.Values(:,1);
        
%% 3) Time strings  
    if Time_SW(1) ~= Time_HD(1) | Time_SW(end) ~= Time_HD(end)
       error(' Start/end time of SW/HD are not the same.')
    end
        
    [y0,m0,d0] = datevec(Time_SW(1));
    [y1,m1,d1] = datevec(Time_SW(end));        
        
%% 4) Set inputs to m_structure
    xyz        = [Lon Lat Depth];
    dt_SW      = (Time_SW(2)-Time_SW(1))*24*60;
    ttt_SW     = [datenum([y0 m0 d0]) datenum([y1 m1 d1]) dt_SW dt_SW];
    dt_HD      = (Time_HD(2)-Time_HD(1))*24*60;
    ttt_HD     = [datenum([y0 m0 d0]) datenum([y1 m1 d1]) dt_HD dt_HD];
        
    % HEADS UP!! HD is every 30 mins, whereas SW is every 60 min ->
    % remove non-matching HD timesteps
    [C,IA,IB]  = intersect(Time_SW,Time_HD);
    Time_HD    = Time_HD(IB);
    WL         = WL(IB);
        
    % Convert to m structure
    Name       = Names{1};
    M_Hm0      = m_structure(Name,xyz,ttt_SW,legend_SW,[Time_SW,Hm0],'Hm0',1,bins_Hm0);
    M_Tp       = m_structure(Name,xyz,ttt_SW,legend_SW,[Time_SW,Tp ],'Tp' ,1,bins_Tp);
    M_T01      = m_structure(Name,xyz,ttt_SW,legend_SW,[Time_SW,T01],'T01',1,bins_T01);
    M_T02      = m_structure(Name,xyz,ttt_SW,legend_SW,[Time_SW,T02],'T02',1,bins_T02);
    M_PWD      = m_structure(Name,xyz,ttt_SW,legend_SW,[Time_SW,PWD],'PWD',1,bins_PWD);
    M_MWD      = m_structure(Name,xyz,ttt_SW,legend_SW,[Time_SW,MWD],'MWD',1,bins_MWD);
    M_WL       = m_structure(Name,xyz,ttt_HD,legend_HD,[Time_HD,WL ],'WL' ,1,bins_WL );
            
%% 5) Storm modes

    % Root dir
    root_dir = pwd;

    % Move to results folder
    cd(out_pth)
               
    % 5.1) Compute storm modes - Forristall
         disp(' ')
         disp('Computing Forristall_H...')
         m_mode(m_directional(M_Hm0,M_MWD), M_T02, M_Tp , M_WL,'IET', 36, 'IEL', 0.7, 'shorttermdist', 'Forristall_H');
         disp('Finished computing Forristall_H')
 
         disp('Computing Forristall_C...')
         m_mode(m_directional(M_Hm0,M_MWD), M_T02, M_T01, M_WL,'IET', 36, 'IEL', 0.7, 'shorttermdist', 'Forristall_C');
         disp('Finished computing Forristall_C')
 
         disp('Computing Forristall_C...')        
         m_mode(m_directional(M_Hm0,M_MWD), M_T02, M_T01,      'IET', 36, 'IEL', 0.7, 'shorttermdist', 'Forristall_C');        
         disp('Finished computing Forristall_C')
 
    % 5.2) Compute storm modes - Glukhovskiy
         disp(' ')
         disp('Computing Glukhovskiy_H...')
         m_mode(m_directional(M_Hm0,M_MWD), M_T02, M_Tp , M_WL,'IET', 36, 'IEL', 0.7, 'shorttermdist', 'Glukhovskiy_H');
         disp('Finished computing Glukhovskiy_H')

    % Go back to root dir
    cd(root_dir)             
