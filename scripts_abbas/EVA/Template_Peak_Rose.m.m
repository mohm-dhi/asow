%"""
%###############################################################################,
%# Created Date:    '2022-12-15'                                                #,
%# Author(s):       Seyed Abbas Jazaeri - s.a.jazaeri@uottawa.ca                #,
%# Encoding:        UTF-8                                                       #,
%# Language:        MATLAB 2022                                                 #,
%# ---------------------------------------------------------------------------- #,
%# This script reads the configuration file and prepares the variables          #,
%#                                                                              #,
%#                                                                              #,
%# ---------------------------------------------------------------------------- #,
%# Copyright (c) (2022) DHI Water & Environment, Inc.                           #,
%################################################################################
%"""

close all
clearvars -except dir variables items labels counter dirs Dists AAPs Estis main_dir EVA_tables
clc

NET.addAssembly('DHI.Mike.Install');
import DHI.Mike.Install.*;
DHI.Mike.Install.MikeImport.SetupLatest({DHI.Mike.Install.MikeProducts.MikeCore});

currd = pwd;

directory = 'D:\Personal\Works\DHI\1._ASOW\Data\03_PROCESSED_MOOD_data\02-OutputLocations\Structs\';

for i = 1:7
    Spd = load(strcat(directory, 'ASOW', string(i), '_all_structs.mat')).asow_params.(variables(counter, 1));
    Dir = load(strcat(directory, 'ASOW', string(i), '_all_structs.mat')).asow_params.(variables(counter, 2));

    xyz = Spd.xyz;
    Spd.label    = convertStringsToChars(labels(counter, 1));
    if (counter == 3 || counter == 9)
        Spd.label = convertStringsToChars(strcat(labels(counter, 1), num2str(xyz(3)*0.25, '%.1f'), "}"));
    elseif (counter == 4 || counter == 10)
        Spd.label = convertStringsToChars(strcat(labels(counter, 1), num2str(xyz(3)*0.5, '%.1f'), "}"));
    elseif (counter == 5 || counter == 11)
        Spd.label = convertStringsToChars(strcat(labels(counter, 1), num2str(xyz(3)*0.75, '%.1f'), "}"));
    elseif (counter == 6 || counter == 12)
        Spd.label = convertStringsToChars(strcat(labels(counter, 1), num2str(xyz(3)*1.0, '%.1f'), "}"));
    end

    Spd.item     = convertStringsToChars(items(1));
    Spd.unit     = 'm/s';
    Spd.isdir    = 0;
    Spd.isspec   = 0;
    Spd.xyz_str = convertStringsToChars(strcat("(", num2str(abs(xyz(1)), '%.2f'), "\circW; ", num2str(xyz(2), '%.2f'), "\circN; d=", num2str(xyz(3), '%.1f'), "mMSL)"));

    Dir.label    = convertStringsToChars(labels(counter, 2));
    if (counter == 3 || counter == 9)
        Dir.label = convertStringsToChars(strcat(labels(counter, 2), num2str(xyz(3)*0.25, '%.1f'), "}"));
    elseif (counter == 4 || counter == 10)
        Dir.label = convertStringsToChars(strcat(labels(counter, 2), num2str(xyz(3)*0.5, '%.1f'), "}"));
    elseif (counter == 5 || counter == 11)
        Dir.label = convertStringsToChars(strcat(labels(counter, 2), num2str(xyz(3)*0.75, '%.1f'), "}"));
    elseif (counter == 6 || counter == 12)
        Dir.label = convertStringsToChars(strcat(labels(counter, 2), num2str(xyz(3)*1.0, '%.1f'), "}"));
    end

    Dir.item     = convertStringsToChars(items(2));
    Dir.unit     = '\circN-to';
    Dir.isdir    = 1;
    Dir.isspec   = 0;
    Dir.xyz_str = convertStringsToChars(strcat("(", num2str(abs(xyz(1)), '%.2f'), "\circW; ", num2str(xyz(2), '%.2f'), "\circN; d=", num2str(xyz(3), '%.1f'), "mMSL)"));
%     Dir.bins     = 0:30:330;


    Spd.T = [1 5 10 25 50 100 500 1000];
    Dir.T = [1 5 10 25 50 100 500 1000];


    mkdir(fullfile(currd,strcat('ASOW', string(i))));
    cd(fullfile(currd,strcat('ASOW', string(i))));

    opt_EVA.EVdist           = Dists(counter, i);
    opt_EVA.EVtype           = 'AAP';
    opt_EVA.EVcrit           = AAPs(counter, i);
    opt_EVA.N_BootStrap      = 100;
    opt_EVA.estimationmethod = convertStringsToChars(Estis(counter, i));
    opt_EVA.plotflag         = [2];
    opt_EVA.EVadju           = 0;
    opt_EVA.ndec             = 1;
    opt_EVA.constfac         = 1e2;

    opt_EVA.optstr           = 'directional';
    m_extreme(m_subseries(Spd,'directional',Dir), opt_EVA);

    cd(currd);
end

clc;