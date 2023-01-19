%%
clear all; clc

addpath(genpath('C:\DHI\Projects\AtlanticShores\Scripts\asow\scripts_juan\normal_analysis'));

%%

spec_dir = 'C:\DHI\Projects\AtlanticShores\Data\2DSpectra\gamma\';
spec_file = 'ASOW_Spectra_ED1f_1979-01-15_0000-2021-12-31_2300_ParamatersTable.mat';

%% convert gamma to m_structure

param_table = load([spec_dir spec_file]);
gamma_vals = param_table.JONSWAP_Gamma;

