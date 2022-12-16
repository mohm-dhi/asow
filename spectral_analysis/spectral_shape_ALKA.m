
clear all; close all; clc

load('\\USDEN1-STOR.DHI.DK\Projects\41806529\02_RAW_MOOD_data\2D Spectra data\Spec0.1_1979-2021_-73.9W39.3N.mat','X')

X.dfs_unit = X.unit;
X.dfs_info.UnitAbbr = 'm^2*s/deg';

%%
X.time = X.time(1:10);
X.data = X.data(1:10,:,:);

IntParam = m_2DspecInt(X);
% IntParam = [Hm0,Tp,T01,T02,Tm10,PWD,MWD,DSD,DSD_p,Spr,Spr_p,PrWD,Tp2,Spr_p2,PWD2];

save('Spectra_Integral_Parameters.mat','IntParam');

aaa
%%
Hm0_range = [0.5 1];
Tp_range = [4 6];
MWD_range = [157.5 202.5];

Id = find(Hm0.data>=Hm0_range(1) & Hm0.data<=Hm0_range(2)...
    & Tp.data>=Tp_range(1) & Tp.data<=Tp_range(2)...
    & MWD.data>=MWD_range(1) & MWD.data<=MWD_range(2));

Sf_mean = nanmean(X.data(Id,:,:),1);

X.data = ED2f_Part.data(1:2,:,:);
X.time = X.time(1:3);
X.data(1,:,:) = Sf_mean;
X.data(2,:,:) = Sf_mean;
X.data(3,:,:) = Sf_mean;

m_spectra_v2(X,'fit',{'JS','OH','PM'},'name',loc,'legend',legend_SWns)
