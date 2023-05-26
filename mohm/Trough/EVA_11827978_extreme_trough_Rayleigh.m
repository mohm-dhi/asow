% Estimate extreme trough (Rayleigh)
% Project: N-7.2 OWF (German North Sea) 11827978
% 2022-11-23 SJA
% Adopted from Hai Long (\\dkcph1-stor.dhi.dk\projects\11827978\SJA\_tempfolder\TroughCalculations\temp12_trough_test_SJA.m)
%

clear; close all; clc;

% addpath(genpath('D:\SJA\AzureGit\potlab_v2\')); %on vnc11827978-84 and\or vnc11827978-55
% addpath(genpath('D:\SJA\AzureGit\potlab_v2\src\'));
% addpath(genpath('D:\SJA\AzureGit\potlab_v2\res\'));
addpath(genpath('\\dkcph1-stor.dhi.dk\projects\11827978\_COMMON_DATA_\MATLAB\potlab_v2\')); %on vnc11827978-84 and\or vnc11827978-55
addpath(genpath('\\dkcph1-stor.dhi.dk\projects\11827978\_COMMON_DATA_\MATLAB\potlab_v2\res'));
addpath(genpath('\\dkcph1-stor.dhi.dk\projects\11827978\_COMMON_DATA_\MATLAB\potlab_v2\src'));
addpath(genpath('\\dkcph1-stor.dhi.dk\projects\11827978\AMSA\'));
addpath(genpath('\\dkcph1-stor.dhi.dk\projects\11827978\SJA\_tempfolder\TroughCalculations\Scripts\'));

% Legends
legend_SW   =  'SW_{DWF2020}';
legend_HD   =  'HD_{DWF2020}';
legend_Wind =  'CFSR';

%% Network Folder
network_folder = '\\DKCPH1-STOR.DHI.DK\Projects\';
project_number = '11827978';
folder_input_files = fullfile(network_folder,project_number,'_COMMON_DATA_','_MODEL_DATA_EXTRACTION_'); %
location_SW = {'37387'};
HD_loc_file = fullfile(folder_input_files,'HD','HD_ElementID.csv');
SW_loc_file = fullfile(folder_input_files,'SW','SW_ElementID.csv'); %'\\dkcph1-stor2.dhi.dk\Projects\11827978_WS\_COMMON_\recover_dwf2020\mood_v1\OtherFiles\SWDWF_ElementID.csv';
[HD_elements,SW_elements,location_HD,ind_HD,ind_SW] = Get_HD_Element(HD_loc_file,SW_loc_file,char(location_SW));
location_HD = {num2str(location_HD)};
close all

% Data Period of the hindcast model
start_date = [1979,7,1,0,0,0];
end_date = [2021,06,30,23,0,0];

% location info - specific to a site and needs to be changed for each site.
WTG_analysis_location = {'P1'};
site_acronym = {'N72'};
hub_heights = [146 165]; %mMSL
CFSR_data_height = [10, 50, 100, 150, 200, hub_heights]; % by default CFSR winds are at 10 m. So make sure that the wind files are also created for the hub height.
xyz_SW = [SW_elements.Centroid_Lat(ind_SW),SW_elements.Centroid_Lon(ind_SW),SW_elements.Water_Depth(ind_SW)*-1];

tt        = [datenum(start_date) datenum(end_date)]; %exclude warmup period at start of HD (and SW)

%% Initialize %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for loc = 1 %[5 1:4 6:8] %1:length(name_all)
    %     point_ID = points(i);
    point_ID = loc;
    name     = WTG_analysis_location{loc};
    xyzd     = xyz_SW(loc,:);
    xyzh     = [xyz_SW(loc,1:2) 10];
    SW_ele = location_SW{loc};
    HD_ele = location_HD{loc};
    
    % files names    
    fname_SW  = fullfile(folder_input_files,'SW',['SW_',SW_ele,'.mat']);
    fname_HD  = fullfile(folder_input_files,'HD',['HD_',HD_ele,'.mat']);
    fname_WS  = fullfile(folder_input_files,'WS_WRA_Aligned',['WS_',SW_ele,'.mat']);
    
    SW = load(fname_SW); HD = load(fname_HD); WS = load(fname_WS);
    %(name,xyz,ttt,legend,fname,item,icol,bins)
    Hm0 = m_structure(name,xyzd,[tt 60 60],legend_SW,[SW.time SW.data],'Hm0',1,0:1:18);
    Tp  = m_structure(name,xyzd,[tt 60 60],legend_SW,[SW.time SW.data],'Tp',2,0:2:26);
    T02 = m_structure(name,xyzd,[tt 60 60],legend_SW,[SW.time SW.data],'T02',4,0:1:18);
    WL  = m_structure(name,xyzd,[tt 30 30],legend_HD,[HD.time HD.data],'WL',1,-1:0.1:1);
    WL  = m_resample(WL,'dt_int',60);

    IET = 48;
    IEL = 0.5;
    outfolder = fullfile(network_folder,project_number,'SJA','04_Analyses','EVA','Trough');
    if ~exist(outfolder,'dir'), mkdir(outfolder); end; cd(outfolder);
    
    % load iPeak from JEVA
    load(fullfile(network_folder,project_number,'AMSA','Struct',sprintf('N72_WTG_%s_10m_wind_H_iPeaks.mat',name)))
    modedata = m_mode_vFLD(Hm0,Tp,T02,WL,...
            'OutputFolder',[outfolder '\modes\'],'plotflag',0,'output_types',1:2,...
            'iPeaks',iPeaks,'intereventtime',IET,'intereventlevel',IEL,...
            'shorttermdist','Rayleigh_T','min_Hm0_peak',0,'return_only_maximum_mode',0);
   
    opt.EVdist              = 4;
    opt.EVcrit              = 1;
    opt.EVtype              = 'AMP';
    opt.N_BootStrap         = 0;
    opt.estimationmethod    = 'LS';
    opt.plotflag            = 1;
    opt.EVadju              = 0;
    opt.T                   = [1 5 10 50 100 500 1000 10000];
    opt.ndec                 = 1;

    m_extreme_vFLD(modedata,opt);
    
    %do the same for height (for comparison)
    modedata = m_mode_vFLD(Hm0,Tp,T02,WL,...
            'OutputFolder',[outfolder '\modes\'],'plotflag',0,'output_types',1:2,...
            'iPeaks',iPeaks,'intereventtime',IET,'intereventlevel',IEL,...
            'shorttermdist','Rayleigh_H','min_Hm0_peak',0,'return_only_maximum_mode',0);
    m_extreme_vFLD(modedata,opt);
end
