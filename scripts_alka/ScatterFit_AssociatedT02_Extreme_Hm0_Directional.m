%"""
%###############################################################################,
%# Created Date:    '2022-01-25'                                                #,
%# Author(s):       ALKA                                                        #,
%# Encoding:        UTF-8                                                       #,
%# Language:        MATLAB                                                      #,
%# ---------------------------------------------------------------------------- #,
%# This script calculates the fit between Tm02 and Hs to estimate associated    #,
%# period to extreme waves.                                                     #,
%# ---------------------------------------------------------------------------- #,
%# Copyright (c) (2022) DHI Water & Environment, Inc.                           #,
%################################################################################
%"""

close all; clear all; clc
addpath(genpath('c:\Users\alka\Documents\MATLAB\potlab_v2-asow\src\m_tools\'))

currd = pwd;

% Define main directory (either network or local)
MainDir = 'c:\Users\alka\Projects\41806529 - ASOW Metocean\';
MStructDir = [MainDir,'03_PROCESSED_MOOD_data\02-OutputLocations\Structs\'];
ExtremesDir = [MainDir,'08_Results\_met_results\03 Extreme\3.4 Waves\3.4.2 Hsig\'];
OutputDir = [MainDir,'08_Results\_met_results\03 Extreme\3.4 Waves\3.4.8  T02 associated with extreme Hs\'];

cor = colororder;

for i = 1:7
    
    % lazy coding style...
    clearvars -except i cor OutputDir ExtremesDir MStructDir MainDir currd
     
    Hm0 = load(strcat(MStructDir, 'ASOW', string(i), '_all_structs.mat')).asow_params.Hm0_Total;
    T02 = load(strcat(MStructDir, 'ASOW', string(i), '_all_structs.mat')).asow_params.T02_Total;
    MWD = load(strcat(MStructDir, 'ASOW', string(i), '_all_structs.mat')).asow_params.MWD_Total;
    
    Hm0.bins = 0:0.5:12;
    T02.bins = 0:2:30;
    
    Hm0.legend = 'SW_{US-EC}';
    T02.legend = 'SW_{US-EC}';
    
    mkdir(fullfile(OutputDir,strcat('ASOW', string(i)),'Directional'));
    cd(fullfile(OutputDir,strcat('ASOW', string(i)),'Directional'));
    
    % calculate fit Hm0 x T02
    quant = [0.05 0.5 0.95];
    fitOut = m_scatter(Hm0,T02,'directional',MWD,'density','plot_std','fit_func','Poly','fit_func_Q','Poly',...
        'Xmin_frac',0.95,...
        'quantiles',quant,'notable','closefig',0);
    
    % Save fit
    fitOut.fit_type = 'Poly';
    fitOut.X = 'Hm0_Total';
    fitOut.Y = 'T02_Total';
    fitOut.ttt = Hm0.ttt;
    fitOut.QuantileP = quant;
    
    save(strcat('ASOW', string(i),'_','Fit_',fitOut.X,'_',fitOut.Y,'.mat'),'fitOut');
    
    % load extremes
    ext_fname = strcat(ExtremesDir,'P',string(i),'\Directional\ASOW',string(i),'_Extreme_Hm0_Data_(1979-01-15_2021-12-31)_AAP_3.00_3_SF=1_LS_Directionaldirectional.txt');
    tab = readtable(ext_fname);
    RP = table2array(tab(1,2:end));
    
    Hm0_ext0 = table2array(tab(2:end,2:end));
    
    dir_order = [0:30:330,nan]; % omni-directional is last
    fig_order = [2:13,1]; % omni-directional is first
    
    % reorder fit output due to divergence between extreme table and fit
    fitOut.All = [fitOut.All(:,2:end) fitOut.All(:,1)];
    fitOut.Quantiles = cat(3,fitOut.Quantiles(:,:,2:end),fitOut.Quantiles(:,:,1));
    
    
    Hm0_ext.Lower = Hm0_ext0(1:3:end,:);
    Hm0_ext.Central = Hm0_ext0(2:3:end,:);
    Hm0_ext.Upper = Hm0_ext0(3:3:end,:);
    labels = table2array(tab(3:3:end,1));
    
    % associated T02 are calculated only for central estimate
    AssociatedT02 = struct();
    %%

    for d = 1:length(dir_order)
        
        if isnan(dir_order(d))
            dir_str = 'Omni';
        else
            dir_str = num2str(dir_order(d));
        end
        
        AssociatedT02.LS(d,:) = fitOut.All(1,d) * Hm0_ext.Central(d,:).^fitOut.All(2,d);
        AssociatedT02.Quant05(d,:) = fitOut.Quantiles(1,1,d) * Hm0_ext.Central(d,:).^fitOut.Quantiles(2,1,d);
        AssociatedT02.Quant95(d,:) = fitOut.Quantiles(1,3,d) * Hm0_ext.Central(d,:).^fitOut.Quantiles(2,3,d);
        
        % add extreme values to figure
        gcf = figure(fig_order(d));
        P = plot(Hm0_ext.Central(d,:),AssociatedT02.LS(d,:),'s','DisplayName','Extremes: H_{m0} Central Est./T_{02} LS Fit','LineWidth',1.2,'Color',cor(3,:));
        plot([Hm0_ext.Central(d,:);Hm0_ext.Central(d,:)],[AssociatedT02.Quant05(d,:);AssociatedT02.Quant95(d,:)],'-',...
            'DisplayName','Extremes: H_{m0} Central Est./T_{02} LS Fit','LineWidth',1.0,...
            'Color',P.Color,'HandleVisibility','off');
        
        % Remove 50% Quantile
        handle_quantiles = findall(gcf,'Linestyle','--');
        delete(handle_quantiles(2));
        
        % Extend Quantile Lines
        mask = Hm0_ext.Central(d,:) > max(handle_quantiles(1).XData);
        handle_quantiles(1).XData = [handle_quantiles(1).XData,Hm0_ext.Central(d,mask)];
        handle_quantiles(1).YData = [handle_quantiles(1).YData,AssociatedT02.Quant05(d,mask)];
        
        mask = Hm0_ext.Central(d,:) > max(handle_quantiles(3).XData);
        handle_quantiles(3).XData = [handle_quantiles(3).XData,Hm0_ext.Central(d,mask)];
        handle_quantiles(3).YData = [handle_quantiles(3).YData,AssociatedT02.Quant95(d,mask)];
        
        % save new figure
        print(strcat('ASOW', string(i),'_ScatterFit_Extremes_Hm0_T02_MWD_',dir_str,'.png'),'-dpng');
    end
    
    % prepare output table
    tab_out = nan(size(Hm0_ext.Central,1)*2+2, length(RP)+1);
    tab_out(2,2:end) = RP;
    
    tab_out(3:2:end,2:end) = Hm0_ext.Central;
    tab_out(4:2:end,2:end) = round(AssociatedT02.LS,1);
    
    tab_out = num2cell(tab_out);
    
    tab_out{1,2} = 'RP';
    
    for d = 1:length(dir_order)
        tab_out{d*2-2+3,1} = ['Hm0 - ', labels{d}];
        tab_out{d*2-2+4,1} = ['Ass. T02 - ', strrep(labels{d},'Central Est.','LS Fit')];
    end
    
    T02_ext = [tab_out(1:2,:); tab_out(4:2:end,:)];
    
    % save results
    save(strcat('ASOW', string(i),'_AssociatedT02.mat'),'AssociatedT02','Hm0_ext');
    
    xlswrite(strcat('ASOW', string(i),'_Extremes_Hm0_T02_MergedTable.xlsx'),tab_out)
    xlswrite(strcat('ASOW', string(i),'_Extremes_Hm0_T02.xlsx'),T02_ext)
    
    cd(currd);
    
    close all
    
end

clc;