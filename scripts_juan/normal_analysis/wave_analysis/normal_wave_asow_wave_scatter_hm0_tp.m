%%
clear all; clc

addpath(genpath('C:\DHI\Projects\AtlanticShores\Scripts\asow\scripts_juan\normal_analysis'));

%%

fdir = 'C:\DHI\Projects\AtlanticShores\Data\TimeSeries\Deliverable\CSVFinal\structs\';
odir_scatter = 'C:\DHI\Projects\AtlanticShores\Analysis\Normal\Waves\Results\Scatter\';
main_dir = 'Hm0_Tp';
no_locs = 7;

%% Hm0/Tp Omni

params1 = {'Hm0_Total','Hm0_Sea','Hm0_Swell'};
params2 = {'Tp_Total','Tp_Sea','Tp_Swell'};
% direc = {'PWD_Total','PWD_Sea','PWD_Swell'};
bins_hm0 = {0:2:8, 0:2:8, 0:2:8};
bins_tp = {0:5:20, 0:5:20, 0:5:20};
leg = {'H_{m0,Total}','H_{m0,Sea}','H_{m0,Swell}'};

for i=1:no_locs

    cd([odir_scatter 'P' num2str(i)]);
    mkdir(main_dir);
    mkdir(main_dir,'Omni');
    cd([main_dir '\Omni']);

    % load structure
    fname = [fdir 'ASOW' num2str(i) '_all_structs.mat'];
    load(fname);

    for p = 1:length(params1)

        mkdir(params1{p});
        cd(params1{p});
        
        hm0_struct = asow_params.(params1{p});
        hm0_struct.data(hm0_struct.data<0)=NaN;
        tp_struct = asow_params.(params2{p});
        tp_struct.data(tp_struct.data<0)=NaN;

        % modify struct items for plotting
        hm0_struct.bins = bins_hm0{p};
        tp_struct.bins = bins_tp{p};
        hm0_struct.legend = 'SW_{US-EC}';
        tp_struct.legend = 'SW_{US-EC}';
        hm0_struct.ttt_str_long = ' (1979-01-15–2021-12-31; \Deltat=1h) ';
        hm0_struct.ttt = [datenum('1979-01-15') datenum('2021-12-31') 60];
        hm0_struct.xyz_str =  ['(' num2str(-1*hm0_struct.xyz(1),'%.3f') 'W; ' num2str(hm0_struct.xyz(2),'%.3f') 'N; ' num2str(hm0_struct.xyz(3),'%.1f') 'm' ')'];

        % do statistics
        [fitvals] = m_scatter(hm0_struct,tp_struct,'density','quantiles',[0.05 0.5 0.95],'fit_func_Q','Poly','nolegend');
        save('fitvals.mat','fitvals');
        cd('..');

    end

end

%% Hm0/Tp Directional

params1 = {'Hm0_Total','Hm0_Sea','Hm0_Swell'};
params2 = {'Tp_Total','Tp_Sea','Tp_Swell'};
direc1 = {'MWD_Total','MWD_Sea','MWD_Swell'};
direc2 = {'PWD_Total','PWD_Sea','PWD_Swell'};
bins_hm0 = {0:2:8, 0:2:8, 0:2:8};
bins_tp = {0:5:20, 0:5:20, 0:5:20};

for i=1:no_locs

    cd([odir_scatter 'P' num2str(i)]);
    mkdir(main_dir);
    mkdir(main_dir,'Directional');
    cd([main_dir '\Directional']);

    % load structure
    fname = [fdir 'ASOW' num2str(i) '_all_structs.mat'];
    load(fname);

    for p = 1:length(params1)

        mkdir(params1{p});
        cd(params1{p});
        
        hm0_struct = asow_params.(params1{p});
        hm0_struct.data(hm0_struct.data<0)=NaN;

        tp_struct = asow_params.(params2{p});
        tp_struct.data(tp_struct.data<0)=NaN;

        mwd_struct = asow_params.(direc1{p});
        mwd_struct.data(mwd_struct.data<0)=NaN;

        pwd_struct = asow_params.(direc2{p});
        pwd_struct.data(pwd_struct.data<0)=NaN;

        % modify struct items for plotting
        hm0_struct.bins = bins_hm0{p};
        tp_struct.bins = bins_tp{p};
        hm0_struct.legend = 'SW_{US-EC}';
        tp_struct.legend = 'SW_{US-EC}';
        hm0_struct.ttt_str_long = ' (1979-01-15–2021-12-31; \Deltat=1h) ';
        hm0_struct.ttt = [datenum('1979-01-15') datenum('2021-12-31') 60];
        hm0_struct.xyz_str =  ['(' num2str(-1*hm0_struct.xyz(1),'%.3f') 'W; ' num2str(hm0_struct.xyz(2),'%.3f') 'N; ' num2str(hm0_struct.xyz(3),'%.1f') 'm' ')'];

        % do statistics
        [fitvals] = m_scatter(hm0_struct,tp_struct,'directional',mwd_struct,'density','quantiles',[0.05 0.5 0.95],'fit_func_Q','Poly','nolegend');
        save('fitvals.mat','fitvals');
        cd('..');

    end

end