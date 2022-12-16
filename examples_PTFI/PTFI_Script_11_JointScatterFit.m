%"""
%###############################################################################,
%# Created Date:    '2022-11-22'                                                #,
%# Author(s):       ALKA                                                        #,
%# Encoding:        UTF-8                                                       #,
%# Language:        MATLAB                                                      #,
%# ---------------------------------------------------------------------------- #,
%# ---------------------------------------------------------------------------- #,
%# Copyright (c) (2022) DHI Water & Environment, Inc.                           #,
%################################################################################
%"""
%--------------------------------------------------------------------------
%01.- Preamble
%--------------------------------------------------------------------------
clc; clear all; close all; fclose all;

addpath(genpath('potlab_v2\src\m_tools\'))
addpath('m_tools_local/')
%run('c:\Users\alka\Documents\MATLAB\wafo-master\wafo\initwafo.m')

%% SETTINGS/OPTIONS

% Define main directory (either network or local)
MainDir = '//usden1-nas3/41806379/';
MStructDir = [MainDir,'06_Results/M_STRUCTURES/'];
OutputDir = [MainDir,'07_Figures/JointOperationalFitQ/'];

% labels
labels = {'total','sea','swell'};
label_name = {'Total Spectrum','Windsea','Swell'};

varnames = {'Sign. Wave Height','Hm0';
    'Wave Period, T02','TZ';
    'Mean Wave Direction','MWD';
    'Dir. Stand. Deviation Mean','DSD';
    'Peak Wave Period','TP';
    'Peak Wave Direction - Robust','PWD'};

% seasons
months = {1:12;
    [12,1,2,3];
    [5,6,7,8];
    [4,9,10,11]};

season_names = {'Annual','Roughest Season','Moderately Rough Season','Mildest Season'};

locations = {'FSRU','Pipeline','MOLF South','MOLF North'};

join_pairs = {
    {'WS','Hm0_total','WD'};
    {'WS','Hm0_swell'};
    {'WS','Hm0_sea','WD'};
    {'WS','WL'};
    {'WS','CS'};
    {'Hm0_total','WL'};
    {'Hm0_total','CS'};
    {'Hm0_total','WS','MWD_total'};
    {'Hm0_sea','WS','MWD_sea'};
    {'Hm0_total','TZ_total','MWD_total'};
    {'Hm0_sea','TZ_sea','MWD_sea'};
    {'Hm0_swell','TZ_swell','MWD_swell'};
    {'Hm0_total','TP_total','MWD_total'};
    {'Hm0_sea','TP_sea','MWD_sea'};
    {'Hm0_swell','TP_swell','MWD_swell'};
    };

%%  LOOP THROUGH LOCATIONS
for g = 1:2%length(locations)
    
    loc = locations{g};
    
    if any(strcmp({'FSRU','Pipeline'},loc))
        loc_wind = 'Offshore';
        load_part = 1;
    elseif any(strcmp({'MOLF South','MOLF North'},loc))
        loc_wind = 'MET-08';
        load_part = 0;
    else
        error('Wind location not assigned')
    end
    
    %% LOOP THROUGH PARTITIONS
    if load_part
        for k = 1:length(labels)
            label = labels{k};
            
            %% LOAD WAVES - OFFSHORE
            for n = 1:size(varnames,1)
                load([MStructDir,loc,'_',varnames{n,2},label,'.mat']);
                eval([varnames{n,2},'_',label, '=', varnames{n,2},';']);
                clear(varnames{n,2})
            end
        end
    else
        %% LOAD WAVES - MOLF
        for n = 1:size(varnames,1)
            load([MStructDir,loc,'_',varnames{n,2},'total.mat']);
            eval([varnames{n,2},'_total =', varnames{n,2},';']);
            clear(varnames{n,2})
        end
    end
    
    %% LOAD WL
    load([MStructDir,loc,'_WL_Combined.mat']);
    WL = WL_Combined;
    clear WL
    
    %% LOAD WIND
    load([MStructDir,loc_wind,'_WS.mat']);
    load([MStructDir,loc_wind,'_WD.mat']);
    
    %% LOAD CURRENTS
    load([MStructDir,loc,'_CS.mat']);
    load([MStructDir,loc,'_CD.mat']);
    
    %% LOOP THROUGH SEASON
    for m = 1:length(season_names)
        Season = season_names{m};
        
        if ~isfolder([OutputDir,loc,'/',Season])
            mkdir([OutputDir,loc,'/',Season])
        end
        cd([OutputDir,loc,'/',Season])
        
        %% LOOP THROUGH PAIRS
        for p = 10:length(join_pairs)
            
            if ~exist(join_pairs{p}{1},'var') || ~exist(join_pairs{p}{2},'var')
                disp('Not found')
                disp(join_pairs{p})
                continue
            end
            
            S = struct();
            S.(join_pairs{p}{1}) = eval([join_pairs{p}{1}]);
            S.(join_pairs{p}{2}) = eval([join_pairs{p}{2}]);
            
            % add directional
            if length(join_pairs{p}) == 3
                S.Dir = eval([join_pairs{p}{3}]);
            end
            
            % ENSURE EVERYTHING IS IN THE SAME TIME
            vars = fieldnames(S);
            len = [];
            for j = 1:length(vars)
                len = [len; length(S.(vars{j}).time)];
            end
            
            while length(unique(len)) > 1
                for j = 1:length(vars)-1
                    
                    X = S.(vars{j});
                    Y = S.(vars{j+1});
                    
                    [X,Y] = m_sync(X,Y,'maxgap',10,'method','nearest');
                    
                    S.(vars{j}) = X;
                    S.(vars{j+1}) = Y;
                    
                end
                
                len = [];
                for j = 1:length(vars)
                    len = [len; length(S.(vars{j}).time)];
                end
            end
            
            for j = 1:length(vars)
                X = S.(vars{j});
                disp(length(X.time))
                S.(vars{j}).name = loc;
            end
            
            %% SPLIT SEASONS
            % get months
            mlist = datevec(S.(vars{j}).time);
            mlist = mlist(:,2);
            
            % mask (data to remove)
            mask = ~ismember(mlist, months{m});
            
            for j = 1:length(vars)
                S.(vars{j}).data(mask) = nan;
            end
            
            % change title
            if m > 1
                S.(vars{j}).ttt_str = strrep(S.(vars{j}).ttt_str,'(',['(', Season,', ']);
            end
            
            if contains(vars{2},'TZ') | contains(vars{2},'TP')
                legloc='SE';
            else
                legloc='NW';
            end
            
            max(S.(vars{2}).data)
            
            if isfield(S,'Dir')
                fitOut = m_scatter_local_PTFI(S.(vars{1}),S.(vars{2}),'density','directional',S.Dir,'plot_std','fit','Poly','fitQ','Poly','quantiles',[0.1 0.9],'dir_str_labels','notable','legloc',legloc);
            else
                fitOut = m_scatter_local_PTFI(S.(vars{1}),S.(vars{2}),'density','plot_std','fit','Poly','fitQ','Poly','quantiles',[0.1 0.9],'notable','legloc',legloc);
            end
            
            % Save fit
            fitOut.fit_type = 'Poly';
            fitOut.X = join_pairs{p}{1};
            fitOut.Y = join_pairs{p}{2};
            fitOut.loc = loc;
            fitOut.ttt = X.ttt;
            save([loc,'_','Fit_',fitOut.X,'_',fitOut.Y,'.mat'],'fitOut');
 
        end
    end
end

