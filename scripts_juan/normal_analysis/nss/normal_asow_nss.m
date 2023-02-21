%%
clear all; clc

addpath(genpath('C:\DHI\Projects\AtlanticShores\Scripts\asow\scripts_juan\normal_analysis'));

%%

spec_dir = 'C:\DHI\Projects\AtlanticShores\Data\2DSpectra\gamma\';
spec_file = 'ASOW_Spectra_ED1f_1979-01-15_0000-2021-12-31_2300_ParamatersTable.mat';

fdir = 'C:\DHI\Projects\AtlanticShores\Data\TimeSeries\Deliverable\CSVFinal\structs\';
no_locs = 7;

%% convert gamma to m_structure

param_table = load([spec_dir spec_file]);
gamma_vals = [datenum(param_table.T.Time(:,:),'dd-mm-yyyy HH:MM') param_table.T.JONSWAP_Gamma];
xyz = [-73.9 39.3 00];
ttt = [datenum([1979 01 15 00 00 00]) datenum([2021 12 31 23 00 00]) 60];

gamma_struct = m_structure('Loc008',xyz,ttt,'Gamma',gamma_vals,'Gamma',1,1:0.25:5);

%%

for i=1:no_locs

    % load structure
    fname = [fdir 'ASOW' num2str(i) '_all_structs.mat'];
    load(fname);

    % make directional subseries for NSS
    WD_hh = asow_params.WindDir_Hub;
    WS_hh = m_subseries(asow_params.WindSpd_Hub,'directional',WD_hh);
    Hm0 = m_subseries(asow_params.Hm0_Total,'directional',WD_hh);
    Tp = m_subseries(asow_params.Tp_Total,'directional',WD_hh);
    T02 = m_subseries(asow_params.T02_Total,'directional',WD_hh);
    Gamma = m_subseries(gamma_struct,'directional',WD_hh);
    WL = m_subseries(asow_params.WL_Total,'directional',WD_hh);
    CS_sfc = m_subseries(asow_params.CurSpd_Surf_Total,'directional',WD_hh);
    CS_btm = m_subseries(asow_params.CurSpd_Bot_Total,'directional',WD_hh);

    % calculate NSS
    [NSS] = m_NSS(WS_hh,WD_hh,Hm0,Tp,T02,Gamma,WL,CS_sfc,CS_btm,3,28,i);

end
