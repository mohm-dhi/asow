%%
% this script reads the concatenated structure files and changes the
% colormap attribute to "colors_DHI"

clear all; clc

addpath(genpath('C:\DHI\Projects\AtlanticShores\Scripts\asow\scripts_juan\'));

%%

struct_dir = 'C:\DHI\Projects\AtlanticShores\Data\TimeSeries\Deliverable\CSVFinal\structs\';
colormap_dir = 'C:\Users\jngz\MATLAB\Projects\potlab_v2\src\m_tools\utilities\colors\';
colormap_name = 'colors_DHI';
no_locs = 7;

%%

% load colormap
load([colormap_dir colormap_name]);

% change attributes on each strcut field
for i = 1:no_locs

    load([struct_dir 'ASOW' num2str(i) '_all_structs.mat']);
    fields = fieldnames(asow_params);

    for f = 1:length(fields)

        asow_params.(fields{f}).ColorMap = ColorMap;
        asow_params.(fields{f}).ColorOrder = ColorOrder;

    end

    save([struct_dir 'ASOW' num2str(i) '_all_structs.mat'],'asow_params','-mat');

end

%% correction for Hm0_Sea = 0

% change attributes on each strcut field
for i = 1:no_locs

    load([struct_dir 'ASOW' num2str(i) '_all_structs.mat']);

    % find Hm0_Sea is zero
    ii = find(asow_params.Hm0_Sea.data==0);

    % assign *_Sea to NaN
    

    save([struct_dir 'ASOW' num2str(i) '_all_structs.mat'],'asow_params','-mat');

end




a = [asow_params.Hm0_Sea.data(ii) asow_params.Tp_Sea.data(ii) asow_params.MWD_Sea.data(ii)];
