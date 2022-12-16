fdir_winds = 'C:\DHI\Projects\AtlanticShores\Data\CFSRWind\';
fdir_wl = 'C:\DHI\Projects\AtlanticShores\Data\TidalAnalysis\wl\MLLW\';
fdir_uv = 'C:\DHI\Projects\AtlanticShores\Data\TidalAnalysis\uv\';
fdir_waves = 'C:\DHI\Projects\AtlanticShores\Data\Points to be delivered SW_HD';
loc_no = 7;

% wind names

wind_lbl = {'WindSpd_10m','WindDir_10m','WindSpd_Hub','WindDir_Hub'};

% water level names

wl_lbl = {'WL_Total','WL_Tide','WL_Res'};

% depth-averaged names

uvd_lbl = {'CurDir_DepAvg_Total','CurSpd_DepAvg_Total',...
    'CurDir_DepAvg_Tide','CurSpd_DepAvg_Tide',...
    'CurDir_DepAvg_Res','CurSpd_DepAvg_Res'};

uvd_data = {'CD','CS','CDTde','CSTde','CDRes','CSRes'};

% depth-levels output field names
uvz_lbl = {'CurDir_Surf_Total',...
    'CurSpd_Surf_Total',...
    'CurDir_Surf_Tide',...
    'CurSpd_Surf_Tide',...
    'CurDir_Surf_Res',...
    'CurSpd_Surf_Res',...
    'CurDir_Bot_Total',...
    'CurSpd_Bot_Total',...
    'CurDir_Bot_Tide',...
    'CurSpd_Bot_Tide',...
    'CurDir_Bot_Res',...
    'CurSpd_Bot_Res',...
    'CurDir_Mid_Total',...
    'CurSpd_Mid_Total',...
    'CurDir_Mid_Tide',...
    'CurSpd_Mid_Tide',...
    'CurDir_Mid_Res',...
    'CurSpd_Mid_Res',...
    'CurDir_25_Total',...
    'CurSpd_25_Total',...
    'CurDir_25_Tide',...
    'CurSpd_25_Tide',...
    'CurDir_25_Res',...
    'CurSpd_25_Res',...
    'CurDir_75_Total',...
    'CurSpd_75_Total',...
    'CurDir_75_Tide',...
    'CurSpd_75_Tide',...
    'CurDir_75_Res',...
    'CurSpd_75_Res'};

% depth-levels directory names (input)
uvz_dir = {'level_1','level_5','level_3','level_4','level_2'};
uvz_data = {'D_total','U_total','D_tide','U_tide','D_res','U_res'};

% wave field names (these structures are processed in this script) 
wave_lbl = {'Hm0_Total'
'Tp_Total'
'T01_Total'
'T02_Total'
'PWD_Total'
'MWD_Total'
'DSD_Total'
'Hm0_Sea'
'Tp_Sea'
'T01_Sea'
'T02_Sea'
'PWD_Sea'
'MWD_Sea'
'Hm0_Swell'
'Tp_Swell'
'T01_Swell'
'T02_Swell'
'PWD_Swell'
'MWD_Swell'};

wave_m_items = {'Hm0_Total'
'Tp_Total'
'T01_Total'
'T02_Total'
'PWD_Total'
'MWD_Total'
'DSD_Total'
'Hm0_Sea'
'Tp_Sea'
'T01_Sea'
'T02_Sea'
'PWD_Sea'
'MWD_Sea'
'Hm0_Swell'
'Tp_Swell'
'T01_Swell'
'T02_Swell'
'PWD_Swell'
'MWD_Swell'};

sw_files = {'P1_US_EastCoast_SW_-73.945_39.651_24.6_3051.4_1979-01-01_2021-12-31_';
    'P2_US_EastCoast_SW_-73.95_39.307_29.9_3458.2_1979-01-01_2021-12-31_';
    'P3_US_EastCoast_SW_-74.044_39.201_24.4_2725.6_1979-01-01_2021-12-31_';
    'P4_US_EastCoast_SW_-74.116_39.161_28.5_2678.6_1979-01-01_2021-12-31_';
    'P5_US_EastCoast_SW_-74.118_39.252_26.8_3106.9_1979-01-01_2021-12-31_';
    'P6_US_EastCoast_SW_-74.111_39.357_23.2_3174_1979-01-01_2021-12-31_';
    'P7_US_EastCoast_SW_-74.212_39.284_21.9_3302_1979-01-01_2021-12-31_'};


%% process water level and currents first

asow_params = struct();

% for all locations
for i = 1:loc_no

    name = ['ASOW' num2str(i)];

    % load winds
    for w = 1:length(wind_lbl)

        % read struct
        fname = [fdir_winds name '_' wind_lbl{w} '.mat'];
        temp = load(fname);
        fields = fieldnames(temp.(wind_lbl{w}));

        % copy struct fields
        for f=1:length(fields)
            asow_params(i).(wind_lbl{w}).(fields{f}) = temp.(wind_lbl{w}).(fields{f});
        end

    end

    % load water levels
    for w = 1:length(wl_lbl)

        % read struct
        fname = [fdir_wl name '_' wl_lbl{w} '.mat'];
        temp = load(fname);
        fields = fieldnames(temp.(wl_lbl{w}));

        % copy struct fields
        for f=1:length(fields)
            asow_params(i).(wl_lbl{w}).(fields{f}) = temp.(wl_lbl{w}).(fields{f});
        end

    end

    % load depth averaged currents
    for w = 1:length(uvd_lbl)

        % read struct
        fname = [fdir_uv name '\' name '_current_analysis.mat'];
        temp = load(fname);
        fields = fieldnames(temp.(uvd_data{w}));

        % copy struct fields
        for f=1:length(fields)
            asow_params(i).(uvd_lbl{w}).(fields{f}) = temp.(uvd_data{w}).(fields{f});
        end

    end

    % load depth-level currents
    % iterate per level, then per variable
    % one main counter to go through final field names

    main_lbl_cnt = 1;

    for dl = 1:length(uvz_dir)

        for df = 1:length(uvz_data)

            % read struct
            fname = [fdir_uv name '\' uvz_dir{dl} '\' uvz_data{df} '.mat'];
            temp = load(fname);
            fields = fieldnames(temp.(uvz_data{df}));

            % copy struct fields
            for f=1:length(fields)
                asow_params(i).(uvz_lbl{main_lbl_cnt}).(fields{f}) = temp.(uvz_data{df}).(fields{f});
            end

            main_lbl_cnt = main_lbl_cnt + 1;

        end

    end

end

%% process waves

for i=1:length(loc_no)