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

ttt_sw            = [datenum([1979 01 01 01 00 00]) datenum([2021 12 31 00 00 00]) 60 180];
ttt_hd            = [datenum([1979 01 01 01 00 00]) datenum([2021 12 31 00 00 00]) 30 30];
ttt_wind          = [datenum([1979 01 01 01 00 00]) datenum([2021 12 31 00 00 00]) 60 120];

% Bins
bins_H         = 0:2:22;
bins_Hmax      = 0:2:30;
bins_T         = 0:2:24;
bins_D         = 0:30:360;
bins_WS        = 0:5:50;
bins_CS        = 0:0.1:1;
bins_WL        = -1.4:0.2:1.4;

E05xyz =	[-72.72	39.97	-57];
E06xyz = 	[-73.43	39.55	-34];

E05xyz_str = ' (72.72W;39.97N;-57 mMSL) ';
E06xyz_str = ' (73.43W;39.55N;-34 mMSL) ';

output_directory = 'C:\Workplace\41806529 - Atlantic Shores\Hmax_Tmax_Cmax\';


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

files_Wind = {
    '\\usden1-stor\Projects\41806529\03_PROCESSED_MOOD_data\02-OutputLocations\CFSR_Wind\ASOW1_WindSpd_Hub.mat'
    '\\usden1-stor\Projects\41806529\03_PROCESSED_MOOD_data\02-OutputLocations\CFSR_Wind\ASOW2_WindSpd_Hub.mat'
    '\\usden1-stor\Projects\41806529\03_PROCESSED_MOOD_data\02-OutputLocations\CFSR_Wind\ASOW3_WindSpd_Hub.mat'
    '\\usden1-stor\Projects\41806529\03_PROCESSED_MOOD_data\02-OutputLocations\CFSR_Wind\ASOW4_WindSpd_Hub.mat'
    '\\usden1-stor\Projects\41806529\03_PROCESSED_MOOD_data\02-OutputLocations\CFSR_Wind\ASOW5_WindSpd_Hub.mat'
    '\\usden1-stor\Projects\41806529\03_PROCESSED_MOOD_data\02-OutputLocations\CFSR_Wind\ASOW6_WindSpd_Hub.mat'
    '\\usden1-stor\Projects\41806529\03_PROCESSED_MOOD_data\02-OutputLocations\CFSR_Wind\ASOW7_WindSpd_Hub.mat'
};

Names               = {'ASOW1', 'ASOW2', 'ASOW3', 'ASOW4', 'ASOW5', 'ASOW6', 'ASOW7'};

Options_Hm0 = {
3 3 'LS'
3 3 'LS'
3 3 'LS'
3 3 'LS'
3 3 'LS'
3 3 'LS'
3 3 'LS'  
};

Options_Wind = {
3 3 'ML'
3 3 'ML'
3 3 'ML'
3 3 'ML'
3 3 'ML'
3 3 'LS'
3 3 'LS'  
};

Options_WL =  {
2 3 'LS'
2 3 'LS'
2 3 'LS'
2 3 'LS'
2 3 'LS'
2 3 'LS'
2 3 'LS'
};

Options_CS = {
3 3 'LS'
3 4 'LS'
3 3 'LS'
3 3 'LS'
3 3 'LS'
3 3 'LS'
3 3 'LS'  
};


for f = 1 : 7 %length(files_SW)
file_SW             = files_SW{f};
file_HD             = files_HD{f};
file_Wind           = files_Wind{f};
Name                = Names{f};


% define bins and site name

[Lon, Lat, Depth]   = get_xyz(file_SW);
xyz                 = [Lon Lat Depth];
legend_SW           = 'SW_{US-EC}';
legend_HD           = 'HD_{US-EC}';

bins_Hm0            = 0:1:12;
bins_CS             = 0:0.1:1;
bins_WS             = 0:5:55;
bins_WL             = -1.6:0.2:1.6;

bins_H               = 0:1:15;
bins_T               = 0:1:27;
bins_D               = 0:30:360;
%% Structures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Hm0 = m_structure(Name , xyz , ttt_sw , legend_SW, file_SW ,'Hm0' , 1 , bins_Hm0);
Hm0.data = naninterp(Hm0.data);

wind_data = load(file_Wind); wind_data = wind_data.WindSpd_Hub;
WindHub = m_structure(Name , xyz , ttt_wind , 'CFSR at Hub Height', [wind_data.time wind_data.data] ,'WS' , 1 , bins_WS);

CS = m_structure(Name , xyz , ttt_hd , legend_HD, file_HD ,'CS' , 2 , bins_CS);
WL = m_structure(Name , xyz , ttt_hd , legend_HD, file_HD ,'WL_{total}' , 1 , bins_WL);
CS = m_resample(CS, 'dt_int', 60);
WL = m_resample(WL, 'dt_int', 60);

%% 2 Options
opt.T = [1,10,50,100, 1000];   
opt.EVdist = 3;        
opt.EVtype = 'AAP';    
opt.EVcrit = 3;        
opt.thetalim = [-80,80];  
opt.F_cond = 0.9;   
opt.g0 = 0;    
opt.quant  = [0.05,0.5,0.95];

WL_opt = opt;
WL_opt.EVdist = Options_WL{f, 1};
WL_opt.EVcrit = Options_WL{f, 1};
WL_opt.estimationmethod = Options_WL{f, 3}; 
WL_opt.thetalim =[-10,10];

CS_opt = opt;
CS_opt.EVdist = Options_CS{f, 1};
CS_opt.EVcrit = Options_CS{f, 1};
CS_opt.estimationmethod = Options_CS{f, 3}; 

Wind_opt = opt;
Wind_opt.EVdist = Options_Wind{f, 1};
Wind_opt.EVcrit = Options_Wind{f, 1};
Wind_opt.estimationmethod = Options_Wind{f, 3};

Hm0_opt = opt;
Hm0_opt.EVdist = Options_Hm0{f, 1};
Hm0_opt.EVcrit = Options_Hm0{f, 1};
Hm0_opt.estimationmethod = Options_Hm0{f, 3};


root_dir = pwd;
Hm0_WL_dir = fullfile(root_dir, 'Hm0_WL');
if ~isfolder(Hm0_WL_dir); mkdir(Hm0_WL_dir); end
cd(Hm0_WL_dir)
if ~isfolder(fullfile(Hm0_WL_dir, ['P' num2str(f)])); mkdir(fullfile(Hm0_WL_dir, ['P' num2str(f)])); end
cd(fullfile(Hm0_WL_dir, ['P' num2str(f)]));
out_strc =  m_conditional_JPA_1plot(Hm0,Hm0_opt,WL,WL_opt);
out_strc =  m_conditional_JPA(Hm0,Hm0_opt,WL,WL_opt);
% m_JPA_plot(out_strc);
csvwrite('Hm0_WL.csv',out_strc{2}.JPA.X_quant);
cd(root_dir)

Hm0_CS_dir = fullfile(root_dir, 'Hm0_CS');
if ~isfolder(Hm0_CS_dir); mkdir(Hm0_CS_dir); end
cd (Hm0_CS_dir)
if ~isfolder(fullfile(Hm0_CS_dir, ['P' num2str(f)])); mkdir(fullfile(Hm0_CS_dir, ['P' num2str(f)])); end
cd(fullfile(Hm0_CS_dir, ['P' num2str(f)]));
out_strc =  m_conditional_JPA_1plot(Hm0,Hm0_opt,CS,CS_opt);
out_strc =  m_conditional_JPA(Hm0,Hm0_opt,CS,CS_opt);
% m_JPA_plot(out_strc);
csvwrite('Hm0_CS.csv',out_strc{2}.JPA.X_quant);
cd(root_dir)


WindHub_Hm0_dir = fullfile(root_dir, 'WindHub_Hm0');
if ~isfolder(WindHub_Hm0_dir); mkdir(WindHub_Hm0_dir); end
cd (WindHub_Hm0_dir)
if ~isfolder(fullfile(WindHub_Hm0_dir, ['P' num2str(f)])); mkdir(fullfile(WindHub_Hm0_dir, ['P' num2str(f)])); end
cd(fullfile(WindHub_Hm0_dir, ['P' num2str(f)]));

out_strc =  m_conditional_JPA_1plot(WindHub,Wind_opt, Hm0,Hm0_opt);
out_strc =  m_conditional_JPA(WindHub,Wind_opt, Hm0,Hm0_opt);
% m_JPA_plot(out_strc);
csvwrite('WindHub_Hm0.csv',out_strc{2}.JPA.X_quant);

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
