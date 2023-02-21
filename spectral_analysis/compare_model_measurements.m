% this script compares wave spectra extracted from MOOD
% with data provided by ASOW/Fugro
%
% ALKA
% Jan 2023

clear all; close all; clc

% load color pallete from m_tools
load('c:\Users\alka\Documents\MATLAB\potlab_v2-asow\src\m_tools\utilities\colors\colors_DHI.mat');

%% PATHS/FILES/SETTINGS
% Main directory with the same folder structure as the datastore
% main_path = '\\USDEN1-STOR.DHI.DK\Projects\41806529\';
main_path = 'c:\Users\alka\Projects\41806529 - ASOW Metocean\';

% Spectra measurements folder
meas_path = [main_path,'01_Client_supplied_data\3 Data\Data\Loc6a\January 2023 Update Spectra Data\Inputs_01122023\Fugro - Spectral Data\DataListing_ASOW\'];

% Spectra files
% the spectra files are for slightly different locations
% check log_file (ASOWBuoy_SpectraDirectionalData_LogSheet.xlsx)
meas_fnames = {'Spectra_SV1'
    'Spectra_SV2_seg1';
    'Spectra_SV2_seg2';
    'Spectra_SV2_seg3';
    'Spectra_SV2_seg4'};

% MOOD file
mood_fname = [main_path,'02_RAW_MOOD_data\2D Spectra data\Spec0.1_1979-2021_-73.9W39.3N_updated_mohm.mat'];

% model bin directional width
% this is important to calculate 1D spectra because model output is in m2/Hz/deg
ddir = 22.5;

% Output Path
fig_path = [main_path,'08_Results\Obs vs Model Comparison\Loc6a-WaveSpectra\'];


%% LOOP PER MEASUREMENT CAMPAING
for n = 1:length(meas_fnames)
    
    %% LOAD DATA - MEASUREMENTS
    % Spectra measurements
    M = load([meas_path,meas_fnames{n}]);
    
    % load and merge all
    % this was an initial idea - superseeded now!
    %     if n == 1
    %         M = load([meas_path,meas_fnames{n}]);
    %     else
    %         M0 = load([meas_path,meas_fnames{n}]);
    %         M.Data = cat(1, M.Data, M0.Data);
    %         M.Time = cat(2, M.Time, M0.Time);
    %     end
    
    % minor clean up
    M.Time = M.Time';
    clear M0
    
    % frequency and direction (from the original LIS files)
    M.freq = 0.04:0.01:0.5;
    M.direc = 30:30:360;
    
    %% LOAD DATA - MODEL
    % MOOD file
    load(mood_fname)
    
    % rename freq/direction for consistency
    X.freq = X.xaxis;
    X.direc = X.yaxis;
    
    %% CROP MODEL DATA TO MEASUREMENT INTERVAL
    mask = (X.time >= min(M.Time)) & (X.time <= max(M.Time));
    X.data = X.data(mask,:,:);
    X.time = X.time(mask);
    
    datestr_start = datestr(min(M.Time),'yyyy/mm/dd');
    datestr_end = datestr(max(M.Time),'yyyy/mm/dd');
    
    %% CALCULATE 1D SPECTRA
    M.Spec1D = sum(M.Data,3);
    X.Spec1D = sum(X.data,3) * ddir; % need to multiply by directional bin with
    
    % averaged over time
    M.Spec1D_Avg = mean(M.Spec1D,1);
    X.Spec1D_Avg = mean(X.Spec1D,1);
    
    %% PLOT AVERAGE SPECTRA
    figure('Color','w','Units','Inches','Position',[1 1 8 3])
    tiledlayout(1,1,'Padding','tight','TileSpacing','tight');
    nexttile;
    hold on
    plot(M.freq, M.Spec1D_Avg, 'LineWidth',1.0, 'DisplayName','Measurements','Color',ColorOrder(1,:));
    plot(X.freq, X.Spec1D_Avg, '--', 'LineWidth',1.2, 'DisplayName','Model','Color',ColorOrder(2,:));
    
    % axis properties
    xlabel('Frequency [Hz]');
    ylabel('Spectral Energy [m^2/Hz]');
    title({'Averaged Spectral Energy'; sprintf('%s - %s', datestr_start, datestr_end)});
    grid on; box on;
    legend('location','NE')
    
    print([fig_path,'Average_1D_Spec_',meas_fnames{n},'.png'],'-dpng','-r300');

    %% PLOT AVERAGE SPECTRA - SQUARE
    figure('Color','w','Units','Inches','Position',[1 1 5 5])
    tiledlayout(1,1,'Padding','tight','TileSpacing','tight');
    nexttile;
    hold on
    plot(M.freq, M.Spec1D_Avg, 'LineWidth',1.0, 'DisplayName','Measurements','Color',ColorOrder(1,:));
    plot(X.freq, X.Spec1D_Avg, '--', 'LineWidth',1.2, 'DisplayName','Model','Color',ColorOrder(2,:));
    
    % axis properties
    xlabel('Frequency [Hz]');
    ylabel('Spectral Energy [m^2/Hz]');
    title({'Averaged Spectral Energy'; sprintf('%s - %s', datestr_start, datestr_end)});
    grid on; box on;
    legend('location','NE')
    ylim([0, 1.3])
    
    print([fig_path,'Average_1D_Spec_',meas_fnames{n},'_SquareAspectRatio.png'],'-dpng','-r300');
   
    %% PLOT AVERAGE SPECTRA - 2 Columns
    figure('Color','w','Units','Inches','Position',[1 1 5 4])
    tiledlayout(1,1,'Padding','tight','TileSpacing','tight');
    nexttile;
    hold on
    plot(M.freq, M.Spec1D_Avg, 'LineWidth',1.0, 'DisplayName','Measurements','Color',ColorOrder(1,:));
    plot(X.freq, X.Spec1D_Avg, '--', 'LineWidth',1.2, 'DisplayName','Model','Color',ColorOrder(2,:));
    
    % axis properties
    xlabel('Frequency [Hz]');
    ylabel('Spectral Energy [m^2/Hz]');
    title({'Averaged Spectral Energy'; sprintf('%s - %s', datestr_start, datestr_end)});
    grid on; box on;
    legend('location','NE')
    ylim([0, 1.3])
    
    print([fig_path,'Average_1D_Spec_',meas_fnames{n},'_2ColAspectRatio.png'],'-dpng','-r300');
    
    %% PLOT HOVMOLLER DIAGRAM
    figure('Color','w','Units','Inches','Position',[1 1 8 5]) ;
    tiledlayout(2,1,'Padding','tight','TileSpacing','tight');
    ax = nexttile;
    pcolor(datetime(M.Time,'convertFrom','datenum'), M.freq, M.Spec1D');
    
    ax(2) = nexttile;
    pcolor(datetime(X.time,'convertFrom','datenum'), X.freq, X.Spec1D');
    
    % colormap
    cmap = pink(20);
    cmap = flipud(cmap);
    
    % colour properties
    for k = 1:length(ax)
        shading(ax(k),'flat');
        caxis(ax(k),[0,10])
        colormap(ax(k),cmap);
        
        % colorbar
        cb = colorbar(ax(k));
        cb.Label.String = 'Spectral Energy [m^2/Hz]';
        cb.TickLabels{end} = ['>',cb.TickLabels{end}];
    end
    
    % axis properties
    for k = 1:length(ax)
        ylim(ax(k),[0,0.6]);
        box(ax(k),'on');
        grid(ax(k),'on');
        ylabel(ax(k),'Frequency [Hz]')
        set(ax(k),'Layer','top')
    end
    
    title(ax(1),'Spectral Energy - Measurements');
    title(ax(2),'Spectral Energy - Model');
    linkaxes(ax)
    
    print([fig_path,'TimeSeries_1D_Spec_',meas_fnames{n},'.png'],'-dpng','-r300');
    
    close all
    
end