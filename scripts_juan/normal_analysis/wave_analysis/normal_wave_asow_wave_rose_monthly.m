%%
clear all; clc

addpath(genpath('C:\DHI\Projects\AtlanticShores\Scripts\asow\scripts_juan\normal_analysis'));


%%

fdir = 'C:\DHI\Projects\AtlanticShores\Data\TimeSeries\Deliverable\CSVFinal\structs\';
odir_rm = 'C:\DHI\Projects\AtlanticShores\Analysis\Normal\Waves\RoseMonthly\';
no_locs = 7;


%% Wave yeights

params = {'Hm0_Total','Hm0_Sea','Hm0_Swell'};
direc = {'MWD_Total','MWD_Sea','MWD_Swell'};
bins = {0:3:12, 0:2:10, 0:2:10};
leg = {'H_{m0,Total}','H_{m0,Sea}','H_{m0,Swell}'};

for i=1:no_locs

    cd(odir_rm);

    % load structure
    fname = [fdir 'ASOW' num2str(i) '_all_structs.mat'];
    load(fname);

    for p = 1:length(params)

        mkdir(params{p});
        cd(params{p});
        mkdir(['ASOW' num2str(i)]);
        cd(['ASOW' num2str(i)])
        wave_struct = asow_params.(params{p});
        wave_struct.data(wave_struct.data<0)=NaN;
        mwd_struct = asow_params.(direc{p});
        mwd_struct.data(mwd_struct.data<0)=NaN;

        % modify struct items for plotting
        wave_struct.bins = bins{p};
        %     wave_struct.unit = 'mMLLW';
        wave_struct.legend = 'SW_{US-EC}';
        wave_struct.ttt_str_long = ' (1979-01-15–2021-12-31; \Deltat=1h) ';
        wave_struct.ttt = [datenum('1979-01-15') datenum('2021-12-31') 60];
        wave_struct.xyz_str =  ['(' num2str(-1*wave_struct.xyz(1),'%.3f') 'W; ' num2str(wave_struct.xyz(2),'%.3f') 'N; ' num2str(wave_struct.xyz(3),'%.1f') 'm' ')'];

        % do statistics
        m_scatter(wave_struct,mwd_struct,'monthly','di',[0 0.01 0.03 0.05 0.1:0.5 0.5:1 1:4]);
        cd('../../');

    end

end

% %% Wave Tp
% 
% params = {'Tp_Total','Tp_Sea','Tp_Swell'};
% direc = {'PWD_Total','PWD_Sea','PWD_Swell'};
% leg = {'T_{p,Total}','T_{p,Sea}','T_{p,Swell}'};
% bins = {0:5:30, 0:5:20, 0:5:30};
% 
% for i=1:no_locs
% 
%     cd(odir_rm);
% 
%     % load structure
%     fname = [fdir 'ASOW' num2str(i) '_all_structs.mat'];
%     load(fname);
% 
%     for p = 1:length(params)
% 
%         mkdir(params{p});
%         cd(params{p});
%         mkdir(['ASOW' num2str(i)]);
%         cd(['ASOW' num2str(i)])
%         wave_struct = asow_params.(params{p});
%         wave_struct.data(wave_struct.data<0)=NaN;
%         pwd_struct = asow_params.(direc{p});
%         pwd_struct.data(pwd_struct.data<0)=NaN;
% 
%         % modify struct items for plotting
%         wave_struct.bins = bins{p};
%         %     wave_struct.unit = 'mMLLW';
%         wave_struct.legend = 'SW_{US-EC}';
%         wave_struct.ttt_str_long = ' (1979-01-15–2021-12-31; \Deltat=1h) ';
%         wave_struct.ttt = [datenum('1979-01-15') datenum('2021-12-31') 60];
%         wave_struct.xyz_str =  ['(' num2str(-1*wave_struct.xyz(1),'%.3f') 'W; ' num2str(wave_struct.xyz(2),'%.3f') 'N; ' num2str(wave_struct.xyz(3),'%.1f') 'm' ')'];
% 
%         % do statistics
%         if p==1 || p==3
%             m_scatter(wave_struct,pwd_struct,'monthly','di',[0 0.5 1:5 8 10 12 15 18]);
%         elseif p==2
%             m_scatter(wave_struct,pwd_struct,'monthly','di',[0 0.5 1:5 8 10]);
%         end
%         cd('../../');
% 
%     end
% 
% end
% 
% %% Wave T01
% 
% params = {'T01_Total','T01_Sea','T01_Swell'};
% direc = {'MWD_Total','MWD_Sea','MWD_Swell'};
% leg = {'T_{01,Total}','T_{01,Sea}','T_{01,Swell}'};
% bins = {0:5:25, 0:5:20, 0:5:25};
% 
% for i=1:no_locs
% 
%     cd(odir_ro);
% 
%     % load structure
%     fname = [fdir 'ASOW' num2str(i) '_all_structs.mat'];
%     load(fname);
% 
%     for p = 3%:length(params)
% 
%         mkdir(params{p});
%         cd(params{p});
%         mkdir(['ASOW' num2str(i)]);
%         cd(['ASOW' num2str(i)])
%         wave_struct = asow_params.(params{p});
%         wave_struct.data(wave_struct.data<0)=NaN;
%         mwd_struct = asow_params.(direc{p});
%         mwd_struct.data(mwd_struct.data<0)=NaN;
% 
%         % modify struct items for plotting
%         wave_struct.bins = bins{p};
%         %     wave_struct.unit = 'mMLLW';
%         wave_struct.legend = 'SW_{US-EC}';
%         wave_struct.ttt_str_long = ' (1979-01-15–2021-12-31; \Deltat=1h) ';
%         wave_struct.ttt = [datenum('1979-01-15') datenum('2021-12-31') 60];
%         wave_struct.xyz_str =  ['(' num2str(-1*wave_struct.xyz(1),'%.3f') 'W; ' num2str(wave_struct.xyz(2),'%.3f') 'N; ' num2str(wave_struct.xyz(3),'%.1f') 'm' ')'];
% 
%         % do statistics
%         if p==1 || p==3
%             m_scatter(wave_struct,mwd_struct,'monthly','di',[0 0.5 1:5 8:10 12]);
%         elseif p==2
%             m_scatter(wave_struct,mwd_struct,'monthly','di',[0 0.1 1:6]);
%         elseif p==3
%             m_scatter(wave_struct,mwd_struct,'monthly','di',[0 0.1 1 5 6:12]);
%         end
%         cd('../../');
% 
%     end
% 
% end
% 
% %% Wave T01, monthly
% 
% params = {'T02_Total','T02_Sea','T02_Swell'};
% direc = {'MWD_Total','MWD_Sea','MWD_Swell'};
% leg = {'T_{02,Total}','T_{02,Sea}','T_{02,Swell}'};
% bins = {0:5:20, 0:5:20, 0:5:20};
% 
% for i=1:no_locs
% 
%     cd(odir_ro);
% 
%     % load structure
%     fname = [fdir 'ASOW' num2str(i) '_all_structs.mat'];
%     load(fname);
% 
%     for p = 1:length(params)
% 
%         mkdir(params{p});
%         cd(params{p});
%         mkdir(['ASOW' num2str(i)]);
%         cd(['ASOW' num2str(i)])
%         wave_struct = asow_params.(params{p});
%         wave_struct.data(wave_struct.data<0)=NaN;
%         mwd_struct = asow_params.(direc{p});
%         mwd_struct.data(mwd_struct.data<0)=NaN;
% 
%         % modify struct items for plotting
%         wave_struct.bins = bins{p};
%         %     wave_struct.unit = 'mMLLW';
%         wave_struct.legend = 'SW_{US-EC}';
%         wave_struct.ttt_str_long = ' (1979-01-15–2021-12-31; \Deltat=1h) ';
%         wave_struct.ttt = [datenum('1979-01-15') datenum('2021-12-31') 60];
%         wave_struct.xyz_str =  ['(' num2str(-1*wave_struct.xyz(1),'%.3f') 'W; ' num2str(wave_struct.xyz(2),'%.3f') 'N; ' num2str(wave_struct.xyz(3),'%.1f') 'm' ')'];
% 
%         % do statistics
%         if p==1 || p==3
%             m_scatter(wave_struct,mwd_struct,'monthly','di',[0 0.5 1:5 8:10]);
%         elseif p==2
%             m_scatter(wave_struct,mwd_struct,'monthly','di',[0 0.1 1:6]);
%         elseif p==3
%             m_scatter(wave_struct,mwd_struct,'monthly','di',[0 0.1 1 5 6:12]);
%         end
%         cd('../../');
% 
%     end
% 
% end