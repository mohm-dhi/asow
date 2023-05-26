% Clean
close all; clear; clc;
% Inst DLL
NET.addAssembly('DHI.Mike.Install');
import DHI.Mike.Install.*;
DHI.Mike.Install.MikeImport.SetupLatest({DHI.Mike.Install.MikeProducts.MikeCore});
NET.addAssembly('DHI.Generic.MikeZero.DFS');
import DHI.Generic.MikeZero.DFS.*;
NET.addAssembly('DHI.Generic.MikeZero.DFS');
import DHI.Generic.MikeZero.DFS.*;
import DHI.Generic.MikeZero.DFS.dfs0.*;

% Add paths
project_name = 'Atlantic_Shores_Offshore_Wind';
output_directory = 'C:\Workplace\41806529 - Atlantic Shores\Hmax_Tmax_Cmax\From Sara LWL with Crest\';
addpath(genpath(output_directory));
%% 0) User input

SP_file             = '\\usden1-stor\Projects\41806529\02_RAW_MOOD_data\2D Spectra data\Spec0.1_1979-2021_-73.9W39.3N_updated_mohm.mat';
SP_file             = 'C:\Workplace\41806529 - Atlantic Shores\Hmax_Tmax_Cmax\Spec0.1_1979-2021_-73.9W39.3N_updated_mohm.mat';
Names               = {'ASOW1', 'ASOW2', 'ASOW3', 'ASOW4', 'ASOW5', 'ASOW6', 'ASOW7'};
files_SW = {
    '\\usden1-stor\Projects\41806529\02_RAW_MOOD_data\Points to be delivered SW_HD\P1_US_EastCoast_SW_-73.945_39.651_24.6_3051.4_1979-01-01_2021-12-31_.dfs0'
    '\\usden1-stor\Projects\41806529\02_RAW_MOOD_data\Points to be delivered SW_HD\P2_US_EastCoast_SW_-73.95_39.307_29.9_3458.2_1979-01-01_2021-12-31_.dfs0'
    '\\usden1-stor\Projects\41806529\02_RAW_MOOD_data\Points to be delivered SW_HD\P3_US_EastCoast_SW_-74.044_39.201_24.4_2725.6_1979-01-01_2021-12-31_.dfs0'
    '\\usden1-stor\Projects\41806529\02_RAW_MOOD_data\Points to be delivered SW_HD\P4_US_EastCoast_SW_-74.116_39.161_28.5_2678.6_1979-01-01_2021-12-31_.dfs0'
    '\\usden1-stor\Projects\41806529\02_RAW_MOOD_data\Points to be delivered SW_HD\P5_US_EastCoast_SW_-74.118_39.252_26.8_3106.9_1979-01-01_2021-12-31_.dfs0'
    '\\usden1-stor\Projects\41806529\02_RAW_MOOD_data\Points to be delivered SW_HD\P6_US_EastCoast_SW_-74.111_39.357_23.2_3174_1979-01-01_2021-12-31_.dfs0'
    '\\usden1-stor\Projects\41806529\02_RAW_MOOD_data\Points to be delivered SW_HD\P7_US_EastCoast_SW_-74.212_39.284_21.9_3302_1979-01-01_2021-12-31_.dfs0'
    };

files_HD = {
    '\\usden1-stor\Projects\41806529\02_RAW_MOOD_data\Points to be delivered SW_HD\P1_US_EastCoast_HD_-73.945_39.651_24.6_3051.4_1979-01-01_2021-12-31_.dfs0'
    '\\usden1-stor\Projects\41806529\02_RAW_MOOD_data\Points to be delivered SW_HD\P2_US_EastCoast_HD_-73.95_39.307_29.9_3458.2_1979-01-01_2021-12-31_.dfs0'
    '\\usden1-stor\Projects\41806529\02_RAW_MOOD_data\Points to be delivered SW_HD\P3_US_EastCoast_HD_-74.044_39.201_24.4_2725.6_1979-01-01_2021-12-31_.dfs0'
    '\\usden1-stor\Projects\41806529\02_RAW_MOOD_data\Points to be delivered SW_HD\P4_US_EastCoast_HD_-74.116_39.161_28.5_2678.6_1979-01-01_2021-12-31_.dfs0'
    '\\usden1-stor\Projects\41806529\02_RAW_MOOD_data\Points to be delivered SW_HD\P5_US_EastCoast_HD_-74.118_39.252_26.8_3106.9_1979-01-01_2021-12-31_.dfs0'
    '\\usden1-stor\Projects\41806529\02_RAW_MOOD_data\Points to be delivered SW_HD\P6_US_EastCoast_HD_-74.111_39.357_23.2_3174_1979-01-01_2021-12-31_.dfs0'
    '\\usden1-stor\Projects\41806529\02_RAW_MOOD_data\Points to be delivered SW_HD\P7_US_EastCoast_HD_-74.212_39.284_21.9_3302_1979-01-01_2021-12-31_.dfs0'
    };
% options: EVdist  EVcrit estimationmethod
Options = {
3 3 'LS'
3 3 'LS'
3 3 'LS'
3 3 'LS'
3 3 'LS'
3 3 'LS'
3 3 'LS'  
};

for f = 1 : 7 %length(files_SW)
file_SW             = files_SW{f};
file_HD             = files_HD{f};
Name                = Names{f};

modes_dir           = sprintf('%s%s\\Modes', output_directory, Name); if ~isfolder(modes_dir); mkdir(modes_dir); end
% out_pth_Hmax        = sprintf('%s%s\\Hmax', output_directory,Name);   if ~isfolder(out_pth_Hmax); mkdir(out_pth_Hmax); end
out_pth_Trough        = sprintf('%s%s\\T10000', output_directory,Name);   if ~isfolder(out_pth_Trough); mkdir(out_pth_Trough); end
% out_pth_Tmax        = sprintf('%s%s\\Tmax', output_directory,Name);   if ~isfolder(out_pth_Tmax); mkdir(out_pth_Tmax); end

% define bins and site name

[Lon, Lat, Depth]   = get_xyz(file_SW);
xyz                 = [Lon Lat Depth];
legend_SW           = 'SW_{US-EC}';
legend_HD           = 'HD_{US-EC}';

bins_Hm0            = 0:1:11;
bins_Tp             = 0:1:30;
bins_T01            = 0:1:20;
bins_T02            = 0:1:20;
bins_WL             = -2.0:0.25:2;
bins_PWD            = 0:30:360;
bins_MWD            = 0:30:360;

bins_H               = 0:1:15;
bins_T               = 0:1:27;
bins_D               = 0:30:360;

bins_ED2            = 0:10 :1000;
%% define EVA parameters
EVAOpt.EVdist           = Options{f, 1};                               % 2. Trunc. Weibull 3. 2-p Weibull 4. Gumbel 5. Exponential
EVAOpt.EVtype           = 'AAP';                                       % Annual Average Peak, quite similar to POT method
EVAOpt.EVcrit           = Options{f, 2};                               % Number of events per year
EVAOpt.N_BootStrap      = 0;                                         % Bootstraping to produce 95% confidence limit
EVAOpt.estimationmethod = Options{f, 3};                               % Least square method to fit extremes onto the distribution
EVAOpt.plotflag         = 1;                                           % 2
EVAOpt.EVadju           = 0;
EVAOpt.T                = [1, 5, 10, 25, 50, 100, 500, 1000, 10000];	       % Return period for second run 1, 5, 10, 25, 50, 100, 500, 1000, 10000
EVAOpt.intereventtime   = 36;	                                       % Inter event time between two extreme events (36 hours)
EVAOpt.intereventlevel  = 0.7;
EVAOpt.ndec             = 1;
EVAOpt.ConfLimits       = [0.025 0.975];                               % 97.5% confidence
EVAOpt.ascii            = '.xlsx';
EVAOpt.constfac         = 1;

%% 1) Load files
if ~exist(strrep(file_SW, '.dfs0', '.mat'), 'file')
    tmp_SW  = load(dfs02mat_dotnet(file_SW));
    tmp_SW = tmp_SW.data;
else
    tmp_SW = load(strrep(file_SW, '.dfs0', '.mat'));
    tmp_SW = tmp_SW.data;
end
if ~exist(strrep(file_HD, '.dfs0', '.mat'), 'file')
    tmp_HD  = load(dfs02mat_dotnet(file_HD));
    tmp_HD = tmp_HD.data;
else
    tmp_HD  = load(strrep(file_HD, '.dfs0', '.mat'));
    tmp_HD = tmp_HD.data;
end
%% 2) Get relevant info
Time_SW    = datenum(datestr(tmp_SW.Time));
Hm0        = tmp_SW.Values(:,1);
Tp         = tmp_SW.Values(:,2);
T01        = tmp_SW.Values(:,3);
T02        = tmp_SW.Values(:,4);
PWD        = tmp_SW.Values(:,5);
MWD        = tmp_SW.Values(:,6);

Time_HD    = datenum(datestr(tmp_HD.Time));
WL         = tmp_HD.Values(:,1);


%% 3) Time strings
if ~strcmp(datestr(Time_HD(1)) ,datestr(Time_SW(1))) || ...
        ~strcmp(datestr(Time_HD(end)) ,datestr(Time_SW(end)))
    error(' Start/end time of SW/HD are not the same.')
end

[y0,m0,d0] = datevec(Time_SW(1));
[y1,m1,d1] = datevec(Time_SW(end));

%% 4) Set inputs to m_structure
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
filename_format = '%s_Mode_SW_{US-EC}_(1979-01-01_-_2021-12-31)_IET=36.0h_IEL=0.70_Rayleigh_%s_%s.mat';

% Move to results folder
cd(modes_dir)

% 5.1) Compute storm modes - Forristall - directional
disp(' ')
disp('Computing Rayleigh_T directional...')
if exist(fullfile(modes_dir, sprintf(filename_format, Name, 'T' ,'MSL_Directional')), 'file') ~= 2
    m_mode_vFLD(m_directional(M_Hm0,M_MWD), M_T02, M_Tp , M_WL, ...
        'IET', 36, 'IEL', 0.7, 'shorttermdist', 'Rayleigh_T', ...
        'output_types',1:2, 'min_Hm0_peak',0,'return_only_maximum_mode',0);
    disp('Finished computing Rayleigh_T')
end

cd(root_dir)

%% 1) Get files in modes_dir
file_list_MSL                 = fullfile(modes_dir, sprintf(filename_format, Name, 'T' ,'MSL_Directional'));

%% 2) Compute Trough
% 2.0) Root dir
root_dir = pwd;
% 2.1) Crests w.r.t. MSL
% Move to dir
cd(out_pth_Trough)

    load(file_list_MSL);
    modedata.unit = '[-]m';
    EVA_T_MSL = m_extreme_vFLD(modedata,EVAOpt);
    
    outname    = [project_name '_Rayleigh_T_MSL_Directional.mat'];
    save(fullfile(out_pth_Trough,outname),'EVA_T_MSL')
cd(root_dir)


end

%% FUNCTIONS
function [lon, lat, depth] = get_xyz(file)
    if contains(file, '_SW_')
        t = split(file, '_SW_'); t = t(end);
        xyz = split(t, '_');
        lon = str2double(xyz{1});
        lat = str2double(xyz{2});
        depth = str2double(xyz{3});
    elseif contains(file, '_HD_')
        t = split(file, '_HD_'); t = t(end);
        xyz = split(t, '_');
        lon = str2double(xyz{1});
        lat = str2double(xyz{2});
        depth = str2double(xyz{3});
    else

    end
end

