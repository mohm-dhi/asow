%%

addpath(genpath('C:\Users\jngz\MATLAB\Projects\potlab\src'));
addpath(genpath('C:\Users\jngz\MATLAB\Projects\potlab\res'));

%% read variable info file

fname_var = 'C:\DHI\Projects\AtlanticShores\Data\TimeSeries\Deliverable\CSV1\asow_csv_variables.csv';
var_map = readtable(fname_var,'NumHeaderLines',0);
var_map = table2cell(var_map);

%% define bins

bin_ws = 0:1:50;
bin_wd = 0:30:330;
bin_wl = 0:0.1:3;
bin_cd = 0:30:330;
bin_cs = 0:0.1:1.5;
bin_wh = 0:0.1:9;
bin_wp = 0:1:25;
bin_wad = 0:30:330;

ttt = [datenum([1979 01 15 00 00 00]) datenum([2021 12 31 23 00 00]) 60];

xyh = [-73.945 39.651 24.6;
    -73.95 39.307 29.9;
    -74.044 39.201 24.4;
    -74.116 39.161 28.5;
    -74.118 39.252 26.8;
    -74.111 39.357 23.2;
    -74.212 39.284 21.9];

loc_no = 7;

%% process files

fdir = 'C:\DHI\Projects\AtlanticShores\Data\TimeSeries\Deliverable\CSV1\';

fname = {'ASOW_P1-DHI-MET-2022 Metocean data_Loc001_x.0';
    'ASOW_P2-DHI-MET-2022 Metocean data_Loc002_x.0';
    'ASOW_P3-DHI-MET-2022 Metocean data_Loc003_x.0';
    'ASOW_P4-DHI-MET-2022 Metocean data_Loc004_x.0';
    'ASOW_P5-DHI-MET-2022 Metocean data_Loc005_x.0';
    'ASOW_P6-DHI-MET-2022 Metocean data_Loc006_x.0';
    'ASOW_P7-DHI-MET-2022 Metocean data_Loc007_x.0'};



% for all locations
for i = 1:loc_no

    asow_params = struct();
    name = ['ASOW' num2str(i)];
    xyz = xyh(i,:);
    xyz_strs = {['(' num2str(xyh(i,1)*-1) 'W;' num2str(xyh(i,2)*-1) 'N;' num2str(xyh(i,3)*-1) 'mMSL)']};
    xyz_str = xyz_strs{1};

    % read csv file

    data_csv = readtable([fdir fname{i} '.csv'],'NumHeaderLines',85);

    % for each variable
    for v = 1:length(var_map)

        data_array = [datenum(data_csv{:,1}) data_csv{:,v+1}];
        asow_params.(var_map{v,1}) = m_structure(name,xyz,ttt,'Data',data_array,var_map{v,1},1,eval(var_map{v,3}));

    end

    %asow_params_out = asow_params(i);
    save(['C:\DHI\Projects\AtlanticShores\Data\TimeSeries\Deliverable\CSV1\structs\' name '_all_structs.mat'], 'asow_params');

end