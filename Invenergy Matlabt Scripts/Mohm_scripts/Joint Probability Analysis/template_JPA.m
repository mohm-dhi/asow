clear ;
clc
% load cmp_data.mat

ttt        = [datenum([1979 01 01 00 00 00]) datenum([2020 12 31 00 00 00]) 60];
% Bins
bins_H         = 0:2:22;
bins_Hmax         = 0:2:30;
bins_T         = 0:2:24;
bins_D         = 0:30:360;
bins_WS         = 0:5:50;
bins_CS         = 0:0.1:1;
bins_WL         = -1.4:0.2:1.4;

E05xyz =	[-72.72	39.97	-57];
E06xyz = 	[-73.43	39.55	-34];

E05xyz_str = ' (72.72W;39.97N;-57 mMSL) ';
E06xyz_str = ' (73.43W;39.55N;-34 mMSL) ';

load ('C:\Users\mohm\OneDrive - DHI\Documents\Project Files\41806038 - Invenergy_NY Bight_Metocean\NYBightData\Hmax\E_05_Mode_SW_{US_East_Coast,CFSR}_(1979-01-01_-_2020-12-31)_IET=36.0h_IEL=0.70_Forristall_H_Directional.mat');
E05Hmax = modedata;
load 'C:\Users\mohm\OneDrive - DHI\Documents\Project Files\41806038 - Invenergy_NY Bight_Metocean\NYBightData\Hmax\E_06_Mode_SW_{US_East_Coast,CFSR}_(1979-01-01_-_2020-12-31)_IET=36.0h_IEL=0.70_Forristall_H_Directional.mat'
E06Hmax = modedata;

load 'C:\Users\mohm\OneDrive - DHI\Documents\Project Files\41806038 - Invenergy_NY Bight_Metocean\NYBightData\Cmax\E05\NY_Bight_Forristall_Cmax_MSL.mat'
load 'C:\Users\mohm\OneDrive - DHI\Documents\Project Files\41806038 - Invenergy_NY Bight_Metocean\NYBightData\Cmax\E05\NY_Bight_Forristall_Cmax_SWL.mat'
E05CmaxMSL = EVA_Cmax_MSL;
E05CmaxSWL = EVA_Cmax_SWL;


E05CmaxMSL.data(:,1:12)=[];
E05CmaxMSL.Modes(:,1:12)=[];
E05CmaxMSL = rmfield(E05CmaxMSL, 'dims');
E05CmaxMSL.label = 'C_{max}';
E05CmaxMSL.vref = '';
E05CmaxMSL.legend = 'CmaxMSL';
E05CmaxMSL.bins = 0:15;

E05CmaxSWL.data(:,1:12)=[];
E05CmaxSWL.Modes(:,1:12)=[];
E05CmaxSWL = rmfield(E05CmaxSWL, 'dims');
E05CmaxSWL.label = 'C_{max}';
E05CmaxSWL.vref = '';
E05CmaxSWL.legend = 'CmaxSWL';
E05CmaxSWL.bins = 0:15;

clear EVA_Cmax_MSL EVA_Cmax_SWL

load 'C:\Users\mohm\OneDrive - DHI\Documents\Project Files\41806038 - Invenergy_NY Bight_Metocean\NYBightData\Cmax\E06\NY_Bight_Forristall_Cmax_MSL.mat'
load 'C:\Users\mohm\OneDrive - DHI\Documents\Project Files\41806038 - Invenergy_NY Bight_Metocean\NYBightData\Cmax\E06\NY_Bight_Forristall_Cmax_SWL.mat'
E06CmaxMSL = EVA_Cmax_MSL;
E06CmaxSWL = EVA_Cmax_SWL;

E06CmaxMSL.data(:,1:12)=[];
E06CmaxMSL.Modes(:,1:12)=[];
E06CmaxMSL = rmfield(E06CmaxMSL, 'dims');
E06CmaxMSL.label = 'C_{max}';
E06CmaxMSL.vref = '';
E06CmaxMSL.legend = 'CmaxMSL';
E06CmaxMSL.bins = 0:15;

E06CmaxSWL.data(:,1:12)=[];
E06CmaxSWL.Modes(:,1:12)=[];
E06CmaxSWL = rmfield(E06CmaxSWL, 'dims');
E06CmaxSWL.label = 'C_{max}';
E06CmaxSWL.vref = '';
E06CmaxSWL.legend = 'CmaxSWL';
E06CmaxSWL.bins = 0:15;


clear modedata;
E05Hmax.data(:, 2:end)=[];
E05Hmax.Modes(:,2:end) = [];
E05Hmax.label = 'H_{max}';
E06Hmax.data(:, 2:end)=[];
E06Hmax.Modes(:,2:end) = [];
E06Hmax.label = 'H_{max}';

load('C:\Users\mohm\OneDrive - DHI\Documents\Project Files\41806038 - Invenergy_NY Bight_Metocean\NYBightData\JWO TidalAnalysis files\E05_current_analysis.mat');
load('C:\Users\mohm\OneDrive - DHI\Documents\Project Files\41806038 - Invenergy_NY Bight_Metocean\NYBightData\JWO TidalAnalysis files\E06_current_analysis.mat');
load('C:\Users\mohm\OneDrive - DHI\Documents\Project Files\41806038 - Invenergy_NY Bight_Metocean\NYBightData\JWO TidalAnalysis files\E05_tidal_analysis.mat')
load('C:\Users\mohm\OneDrive - DHI\Documents\Project Files\41806038 - Invenergy_NY Bight_Metocean\NYBightData\JWO TidalAnalysis files\E06_tidal_analysis.mat')

E05Hm0_file = 'C:\Users\mohm\OneDrive - DHI\Documents\Project Files\41806038 - Invenergy_NY Bight_Metocean\NYBightData\Waves E05\US_East_Coast_Wave_Parameters_Integrated_MIKE_21_Spectral_Wave_Model_DHI_-72.71669_39.96928_NoHurricans.dfs0';
E06Hm0_file = 'C:\Users\mohm\OneDrive - DHI\Documents\Project Files\41806038 - Invenergy_NY Bight_Metocean\NYBightData\Waves E06\US_East_Coast_Wave_Parameters_Integrated_MIKE_21_Spectral_Wave_Model_DHI_-73.42889_39.54677.dfs0';

E05WindHub_file = 'C:\Users\mohm\OneDrive - DHI\Documents\Project Files\41806038 - Invenergy_NY Bight_Metocean\NYBightData\CFSR Hub Height\CFSR Hub 165 E05\CFSR_at_hubHeight_165m.dfs0';
E06WindHub_file = 'C:\Users\mohm\OneDrive - DHI\Documents\Project Files\41806038 - Invenergy_NY Bight_Metocean\NYBightData\CFSR Hub Height\CFSR Hub 165 E06\CFSR_at_hubHeight_165m.dfs0';

%% Structures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E05Hm0 = m_structure('E05' , E05xyz , ttt , 'MIKE SW', E05Hm0_file ,'Hm0' , 1 , bins_H);
E06Hm0 = m_structure('E06' , E06xyz , ttt , 'MIKE SW', E06Hm0_file ,'Hm0' , 1 , bins_H);
E05Hm0.data = naninterp(E05Hm0.data);
E06Hm0.data = naninterp(E06Hm0.data);

E05WindHub = m_structure('E05' , E05xyz , ttt , 'CFSR at Hub Height', E05WindHub_file ,'WS' , 1 , bins_WS);
E06WindHub = m_structure('E06' , E06xyz , ttt , 'CFSR at Hub Height', E06WindHub_file ,'WS' , 1 , bins_WS);



E05Hmax.bins = bins_Hmax;
E06Hmax.bins = bins_Hmax;

E05WL_tot = E05_tidal_analysis.WL;
E05WL_tot.item = 'WL_{total}'; E05WL_tot.label = 'WL_{total}'; E05WL_tot.unit = 'mMSL';
E05WL_tot.xyz = E05xyz;
E05WL_tot.xyz_str = E05xyz_str;
E05WL_tot.bins = bins_WL;


E06WL_tot = E06_tidal_analysis.WL;
E06WL_tot.item = 'WL_{total}'; E06WL_tot.label = 'WL_{total}'; E06WL_tot.unit = 'mMSL';
E06WL_tot.xyz = E06xyz;
E06WL_tot.xyz_str = E06xyz_str;
E06WL_tot.bins = bins_WL;

E05WL_res = E05_tidal_analysis.WLRes;
E05WL_res.item = 'WL_{residual}'; E05WL_res.label = 'WL_{residual}'; E05WL_res.unit = 'mMSL';
E05WL_res.xyz = E05xyz;
E05WL_res.xyz_str = E05xyz_str;


E06WL_res = E06_tidal_analysis.WLRes;
E06WL_res.item = 'WL_{residual}'; E06WL_res.label = 'WL_{residual}'; E06WL_res.unit = 'mMSL';
E06WL_res.xyz = E06xyz;
E06WL_res.xyz_str = E06xyz_str;


E05CS = E05_current_analysis.CS;
E05CS.xyz = E05xyz;
E05CS.xyz_str = E05xyz_str;


E06CS = E06_current_analysis.CS;
E06CS.xyz = E06xyz;
E06CS.xyz_str = E06xyz_str;


E05CS_res = E05_current_analysis.CSRes;
E05CS_res.xyz = E05xyz;
E05CS_res.xyz_str = E05xyz_str;


E06CS_res = E06_current_analysis.CSRes;
E06CS_res.xyz = E06xyz;
E06CS_res.xyz_str = E06xyz_str;

%% 2 Options
opt1.T = [1,10,50,100];   
opt1.EVdist = 3;        
opt1.EVtype = 'AAP';    
opt1.EVcrit = 3;        
opt1.thetalim = [-80,80];  
opt1.F_cond = [0.9];   
opt1.g0 = 0;    
opt1.quant  = [0.05,0.5,0.95];
%% 2.1 WL Options
WL_opt=opt1;
WL_opt.EVdist = 5;
WL_opt.EVtype = 'AAP';
WL_opt.EVcrit = 1;
WL_opt.thetalim =[-10,10];
WL_opt.JPA_sign = 0;  % controls whether marginal distribution of conditioned
%% 2.2 CS options
CS_opt = opt1;
CS_opt.estimationmethod = 'LS';
CS_res_opt = opt1;
CS_res_opt.estimationmethod = 'LS';
%% 2.3 Hm0 options
Hm0E05_opt = opt1;
Hm0E06_opt = opt1;
Hm0E06_opt.EVdist=5;
Hm0E06_opt.estimationmethod = 'LS';

%% 2.4 Crst
CmaxE05_opt = Hm0E05_opt;
CmaxE06_opt = Hm0E06_opt;

%% 2.5 Hmax
HmaxE05 = Hm0E05_opt;
HmaxE06 = Hm0E06_opt;


%% 2.6 WindSpd at Hub heigh options
WindHub_opt = opt1;

%% Run JPA on all WL (no separation into low and high water)
close all,
if ~exist(fullfile(pwd,'EE05_test'), 'dir')
mkdir(fullfile(pwd,'EE05_test'));
end
cd(fullfile(pwd,'EE05_test'));
out_strc =  m_conditional_JPA(E05Hm0,Hm0E05_opt,E05WL_tot,WL_opt);
m_JPA_plot(out_strc);
csvwrite('E05_Hm0_WLtotal.csv',out_strc{2}.JPA.X_quant);

% out_strc =  m_conditional_JPA(E05WL_tot,WL_opt,E05Hm0,Hm0E05_opt);
% m_JPA_plot(out_strc);
% csvwrite('E05_WLtotal_Hm0.csv',out_strc{2}.JPA.X_quant);

out_strc =  m_conditional_JPA(E05Hm0,Hm0E05_opt,E05WL_res,WL_opt);
m_JPA_plot(out_strc);
csvwrite('E05_Hm0_WLres.csv',out_strc{2}.JPA.X_quant);

% out_strc =  m_conditional_JPA(E05WL_res,WL_opt,E05Hm0,Hm0E05_opt);
% m_JPA_plot(out_strc);
% csvwrite('E05_WLres_Hm0.csv',out_strc{2}.JPA.X_quant);

out_strc =  m_conditional_JPA(E05Hmax,opt1,E05WL_tot,WL_opt);
out_strc{2}.JPA.X_quant(2:end,2) = out_strc{1}.R_mx;
m_JPA_plot(out_strc);
csvwrite('E05_Hmax_WLtotal.csv',out_strc{2}.JPA.X_quant);

% out_strc =  m_conditional_JPA(E05WL_tot,WL_opt,E05Hmax,opt1);
% out_strc{2}.JPA.X_quant(2:end,2) = out_strc{1}.R_mx;
% m_JPA_plot(out_strc);
% csvwrite('E05_WLtotal_Hmax.csv',out_strc{2}.JPA.X_quant);

out_strc =  m_conditional_JPA(E05Hmax,opt1,E05WL_res,WL_opt);
out_strc{2}.JPA.X_quant(2:end,2) = out_strc{1}.R_mx;
m_JPA_plot(out_strc);
csvwrite('E05_Hmax_WLres.csv',out_strc{2}.JPA.X_quant);

% out_strc =  m_conditional_JPA(E05WL_res,WL_opt,E05Hmax,opt1);
% out_strc{2}.JPA.X_quant(2:end,2) = out_strc{1}.R_mx;
% m_JPA_plot(out_strc);
% csvwrite('E05_WLres_Hmax.csv',out_strc{2}.JPA.X_quant);

out_strc =  m_conditional_JPA(E05CmaxMSL,CmaxE05_opt,E05WL_tot,WL_opt);
out_strc{2}.JPA.X_quant(2:end,2) = out_strc{1}.R_mx;
m_JPA_plot(out_strc);
csvwrite('E05_CmaxMSL_WLtotal.csv',out_strc{2}.JPA.X_quant);

out_strc =  m_conditional_JPA(E05CmaxSWL,CmaxE05_opt,E05WL_tot,WL_opt);
out_strc{2}.JPA.X_quant(2:end,2) = out_strc{1}.R_mx;
m_JPA_plot(out_strc);
csvwrite('E05_CmaxSWL_WLtotal.csv',out_strc{2}.JPA.X_quant);


out_strc =  m_conditional_JPA(E05CmaxMSL,CmaxE05_opt,E05WL_res,WL_opt);
out_strc{2}.JPA.X_quant(2:end,2) = out_strc{1}.R_mx;
m_JPA_plot(out_strc);
csvwrite('E05_CmaxMSL_WLres.csv',out_strc{2}.JPA.X_quant);

out_strc =  m_conditional_JPA(E05CmaxSWL,CmaxE05_opt,E05WL_res,WL_opt);
out_strc{2}.JPA.X_quant(2:end,2) = out_strc{1}.R_mx;
m_JPA_plot(out_strc);
csvwrite('E05_CmaxSWL_WLres.csv',out_strc{2}.JPA.X_quant);

out_strc =  m_conditional_JPA(E05Hm0,Hm0E05_opt, E05CS,CS_opt);
m_JPA_plot(out_strc);
csvwrite('E05_Hm0_CStotal.csv',out_strc{2}.JPA.X_quant);

% out_strc =  m_conditional_JPA(E05CS,CS_opt,E05Hm0,Hm0E05_opt);
% m_JPA_plot(out_strc);
% csvwrite('E05_CStotal_Hm0.csv',out_strc{2}.JPA.X_quant);


out_strc =  m_conditional_JPA(E05WindHub,WindHub_opt,E05Hm0,Hm0E05_opt);
m_JPA_plot(out_strc);
csvwrite('E05_WindHub_Hm0.csv',out_strc{2}.JPA.X_quant);

out_strc =  m_conditional_JPA(E05Hm0,Hm0E05_opt, E05WindHub,WindHub_opt);
m_JPA_plot(out_strc);
csvwrite('E05_Hm0_WindHub.csv',out_strc{2}.JPA.X_quant);

%% plot for side E06
cd .. 
close all,
if ~exist(fullfile(pwd,'EE06'), 'dir')
mkdir(fullfile(pwd,'EE06'));
end
cd(fullfile(pwd,'EE06'));
out_strc =  m_conditional_JPA(E06Hm0,Hm0E06_opt,E06WL_tot,WL_opt);
m_JPA_plot(out_strc);
csvwrite('E06_Hm0_WLtotal.csv',out_strc{2}.JPA.X_quant);

% out_strc =  m_conditional_JPA(E06WL_tot,WL_opt,E06Hm0,Hm0E06_opt);
% m_JPA_plot(out_strc);
% csvwrite('E06_WLtotal_Hm0.csv',out_strc{2}.JPA.X_quant);

out_strc =  m_conditional_JPA(E06Hm0,Hm0E06_opt,E06WL_res,WL_opt);
m_JPA_plot(out_strc);
csvwrite('E06_Hm0_WLres.csv',out_strc{2}.JPA.X_quant);

% out_strc =  m_conditional_JPA(E06WL_res,WL_opt,E06Hm0,Hm0E06_opt);
% m_JPA_plot(out_strc);
% csvwrite('E06_WLres_Hm0.csv',out_strc{2}.JPA.X_quant);

out_strc =  m_conditional_JPA(E06Hmax,opt1,E06WL_tot,WL_opt);
out_strc{2}.JPA.X_quant(2:end,2) = out_strc{1}.R_mx;
m_JPA_plot(out_strc);
csvwrite('E06_Hmax_WLtotal.csv',out_strc{2}.JPA.X_quant);

% out_strc =  m_conditional_JPA(E06WL_tot,WL_opt,E06Hmax,opt1);
% out_strc{2}.JPA.X_quant(2:end,2) = out_strc{1}.R_mx;
% m_JPA_plot(out_strc);
% csvwrite('E06_WLtotal_Hmax.csv',out_strc{2}.JPA.X_quant);

out_strc =  m_conditional_JPA(E06Hmax,opt1,E06WL_res,WL_opt);
out_strc{2}.JPA.X_quant(2:end,2) = out_strc{1}.R_mx;
m_JPA_plot(out_strc);
csvwrite('E06_Hmax_WLres.csv',out_strc{2}.JPA.X_quant);

% out_strc =  m_conditional_JPA(E06WL_res,WL_opt,E06Hmax,opt1);
% out_strc{2}.JPA.X_quant(2:end,2) = out_strc{1}.R_mx;
% m_JPA_plot(out_strc);
% csvwrite('E06_WLres_Hmax.csv',out_strc{2}.JPA.X_quant);

out_strc =  m_conditional_JPA(E06CmaxMSL,CmaxE06_opt,E06WL_tot,WL_opt);
out_strc{2}.JPA.X_quant(2:end,2) = out_strc{1}.R_mx;
m_JPA_plot(out_strc);
csvwrite('E06_CmaxMSL_WLtotal.csv',out_strc{2}.JPA.X_quant);

out_strc =  m_conditional_JPA(E06CmaxSWL,CmaxE06_opt,E06WL_tot,WL_opt);
out_strc{2}.JPA.X_quant(2:end,2) = out_strc{1}.R_mx;
m_JPA_plot(out_strc);
csvwrite('E06_CmaxSWL_WLtotal.csv',out_strc{2}.JPA.X_quant);


out_strc =  m_conditional_JPA(E06CmaxMSL,CmaxE06_opt,E06WL_res,WL_opt);
out_strc{2}.JPA.X_quant(2:end,2) = out_strc{1}.R_mx;
m_JPA_plot(out_strc);
csvwrite('E06_CmaxMSL_WLres.csv',out_strc{2}.JPA.X_quant);

out_strc =  m_conditional_JPA(E06CmaxSWL,CmaxE06_opt,E06WL_res,WL_opt);
out_strc{2}.JPA.X_quant(2:end,2) = out_strc{1}.R_mx;
m_JPA_plot(out_strc);
csvwrite('E06_CmaxSWL_WLres.csv',out_strc{2}.JPA.X_quant);

out_strc =  m_conditional_JPA(E06Hm0,Hm0E06_opt, E06CS,CS_opt);
m_JPA_plot(out_strc);
csvwrite('E06_Hm0_CStotal.csv',out_strc{2}.JPA.X_quant);

% out_strc =  m_conditional_JPA(E06CS,CS_opt,E06Hm0,Hm0E06_opt);
% m_JPA_plot(out_strc);
% csvwrite('E06_CStotal_Hm0.csv',out_strc{2}.JPA.X_quant);


out_strc =  m_conditional_JPA(E06WindHub,WindHub_opt,E06Hm0,Hm0E06_opt);
m_JPA_plot(out_strc);
csvwrite('E06_WindHub_Hm0.csv',out_strc{2}.JPA.X_quant);

out_strc =  m_conditional_JPA(E06Hm0,Hm0E06_opt, E06WindHub,WindHub_opt);
m_JPA_plot(out_strc);
csvwrite('E06_Hm0_WindHub.csv',out_strc{2}.JPA.X_quant);
