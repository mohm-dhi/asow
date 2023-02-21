%%
% this script fixes Hm0,MWD,TP when "no sea"

clear all; clc

addpath(genpath('C:\DHI\Projects\AtlanticShores\Scripts\asow\scripts_juan\'));

%%

struct_dir = 'C:\DHI\Projects\AtlanticShores\Data\TimeSeries\Deliverable\CSVFinal\structs\';
no_locs = 7;

%% correction for Hm0_Sea = 0

% change attributes on each strcut field
for i = 1:no_locs

    load([struct_dir 'ASOW' num2str(i) '_all_structs.mat']);

    % find Hm0_Sea is zero
    ii = find(asow_params.Hm0_Sea.data==0);

    % assign *_Sea to NaN
    asow_params.Hm0_Sea.data(ii) = NaN;
    asow_params.Tp_Sea.data(ii) = NaN;
    asow_params.MWD_Sea.data(ii) = NaN;

    save([struct_dir 'ASOW' num2str(i) '_all_structs.mat'],'asow_params','-mat');

end

%% modify struct to output Excel spreadsheet

% change attributes on each struct field
for i = 1:no_locs

    load([struct_dir 'ASOW' num2str(i) '_all_structs.mat']);
    fields = fieldnames(asow_params);

    for f = 1:length(fields)

        asow_params.(fields{f}).ascii = '.xlsx';
        asow_params.(fields{f}).name = ['P' num2str(i)];

    end

    save([struct_dir 'ASOW' num2str(i) '_all_structs.mat'],'asow_params','-mat');

end