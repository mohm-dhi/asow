%%

addpath(genpath('\\USDEN1-STOR.DHI.DK\Projects\41806529\_metocean_scripts\potlab_v2\src'));
addpath(genpath('\\USDEN1-STOR.DHI.DK\Projects\41806529\_metocean_scripts\potlab_v2\res'));

% main directory
fdir_tide = 'C:\Users\jngz\OneDrive - DHI\2022\41806529-AtlanticShores\Scripts\FileIO\uv\';
% fdir_cfsr = 'C:\Users\jngz\OneDrive - DHI\2022\41806529-AtlanticShores\Data\CFSR MOOD\OutputLocations\';

fdir_cfsr = 'C:\Users\jngz\OneDrive - DHI\2022\41806529-AtlanticShores\Data\CFSR interpolated\concatenated\';

% cfsr_names = {'P1_Global_Wind_CFSR_286.055_39.651_10_22716.1_1979-01-01_2022-10-01_.dfs0';
%     'P2_Global_Wind_CFSR_286.05_39.307_10_22716.1_1979-01-01_2022-10-01_.dfs0';
%     'P3_Global_Wind_CFSR_285.956_39.201_10_22716.1_1979-01-01_2022-10-01_.dfs0';
%     'P4_Global_Wind_CFSR_285.884_39.161_10_22716.1_1979-01-01_2022-10-01_.dfs0';
%     'P5_Global_Wind_CFSR_285.882_39.252_10_22716.1_1979-01-01_2022-10-01_.dfs0';
%     'P6_Global_Wind_CFSR_285.889_39.357_10_22716.1_1979-01-01_2022-10-01_.dfs0';
%     'P7_Global_Wind_CFSR_285.788_39.284_10_22716.1_1979-01-01_2022-10-01_.dfs0'};

cfsr_names = {'P1_CFSR-global-_1979-01-01_-_2022-12-31.dfs0';
    'P2_CFSR-global-_1979-01-01_-_2022-12-31.dfs0';
    'P3_CFSR-global-_1979-01-01_-_2022-12-31.dfs0';
    'P4_CFSR-global-_1979-01-01_-_2022-12-31.dfs0';
    'P5_CFSR-global-_1979-01-01_-_2022-12-31.dfs0';
    'P6_CFSR-global-_1979-01-01_-_2022-12-31.dfs0';
    'P7_CFSR-global-_1979-01-01_-_2022-12-31.dfs0';};

% lon lat depth
xyh = [-73.945 39.651 24.6;
    -73.95 39.307 29.9;
    -74.044 39.201 24.4;
    -74.116 39.161 28.5;
    -74.118 39.252 26.8;
    -74.111 39.357 23.2;
    -74.212 39.284 21.9];

% depths
site_depth = xyh(:,3);

% bins and dates
bins_D         = 0:30:360; % wind dir
bins_WS         = 0:3:33;  % wind speed
bins_CS         = 0:0.1:1; % current speed
bins_U = -30:0.1:30;
bins_V = -30:0.1:30;
ttt = [datenum([1979 01 01 00 00 00]) datenum([2021 12 31 23 00 00]) 60];

% structs
result = struct;
result_residual = struct;

%%

% for each location
for i=1:length(xyh)

    % descr
    name = ['ASOW' num2str(i)];
    result.(name) = [];
    result_residual.(name) = [];
    depth = site_depth(i);
    xyz = xyh(i,:);
    xyz_strs = {['(' num2str(xyh(i,1)*-1) 'W;' num2str(xyh(i,2)*-1) 'N;' num2str(xyh(i,3)*-1) 'mMSL)']};
    xyz_str = xyz_strs{1};

    % enter tide dir
    cd([fdir_tide name])

    % load tide and cfsr files
    tides = load([name '_current_analysis.mat']);
    cfsr_path = [fdir_cfsr cfsr_names{i}];

    CSTde = tides.CSTde;
    CDTde = tides.CDTde;
    CFSR_U = m_structure(name , xyz , ttt , 'Data' , cfsr_path, 'u-wind' , 1, bins_U);
    CFSR_V = m_structure(name , xyz , ttt , 'Data' , cfsr_path, 'v-wind' , 2, bins_V);
    [WindSpeed,WindDir] = uv2spddir(CFSR_U, CFSR_V, 'invert');
    CFSR_Speed = m_structure(name , xyz , ttt , 'Data' , [CFSR_U.time WindSpeed.data], 'wind speed' , 1, bins_WS);
    CFSR_Dir = m_structure(name , xyz , ttt , 'Data' , [CFSR_U.time WindDir.data], 'wind dir' , 1, bins_D);


    % calculate offsets for desired depths
    fr = [1,0.75,0.5,0.25,0.05];
    z = (depth.*fr) - depth;

    % level 1 for surface
    % level 5 for bottom

    % depth levels
    for zi = 1:length(z)

        % make directory for level
        sdirn = ['level_' num2str(zi)];
        mkdir(sdirn);

        % zb = 0 at bottom
        % zb = depth at surface
        zb = depth + z(zi);

        if zb>0

            tide = timetable(CSTde.datetime, CSTde.data, CDTde.data);
            tide.Properties.VariableNames = {'spd', 'dir'};
            cfsr = timetable(CFSR_Speed.datetime, CFSR_Speed.data, CFSR_Dir.data);
            cfsr.Properties.VariableNames = {'spd', 'dir'};
            tide = retime(tide, cfsr.Time, 'linear');

            % do calculations
            u_tide = tide.spd*(zb/(0.32*depth))^(1/7);
            u_surf = 0.01 * cfsr.spd;

            % make m_structure for tidal current at depth
            U_tide = m_structure(name, xyz, ttt,'CS_{Tide}',[datenum(tide.Time) u_tide tide.dir], 'Current_speed_tide', 1, bins_CS);
            D_tide = m_structure(name, xyz, ttt,'CD_{Tide}',[datenum(tide.Time) u_tide tide.dir], 'Current_direction_tide', 2, bins_D);

            save([sdirn '\U_tide.mat'],'U_tide');
            save([sdirn '\D_tide.mat'],'D_tide');

            % make m_structure for surface current
            U_surf = m_structure(name, xyz, ttt,'CS_{Surface}',[datenum(tide.Time) u_surf cfsr.dir], 'Current_speed_surface', 1, bins_CS);
            D_surf = m_structure(name, xyz, ttt,'CD_{Surface}',[datenum(tide.Time) u_surf cfsr.dir], 'Current_direction_surface', 2, bins_D);

            save([sdirn '\U_surf.mat'],'U_surf');
            save([sdirn '\D_surf.mat'],'D_surf');

            % calculate residual current
            if z(zi)>= -20
                u_res = u_surf .* (1+z(zi)/20);
                % try this too, z has to be positive
                % u_res = u_surf .* (1+z/(depth/2));
            else
                u_res = u_surf * 0;
            end

            % make structures for residual current
            U_res = m_structure(name, xyz, ttt, ['CS_{Residual} z=' num2str(z(zi))], [datenum(tide.Time) u_res], 'Current_speed_residual', 1, bins_CS);
            D_res = m_structure(name, xyz, ttt, ['CD_{Residual} z=' num2str(z(zi))], [datenum(tide.Time) cfsr.dir], 'Current_direction_residual', 1, bins_D);
            U_res.xyz_str = xyz_str;
            D_res.xyz_str = xyz_str;

            % save variables for residuals
            save([sdirn '\U_res.mat'],'U_res');
            save([sdirn '\D_res.mat'],'D_res');

            % split into u,v components
            [u_tide_component, v_tide_component] = spddir2uv(u_tide, tide.dir);
            [u_surf_component, v_surf_component] = spddir2uv(u_res, cfsr.dir, 'invert');

            % calculate total current
            u_tot = u_tide_component + u_surf_component;
            v_tot = v_tide_component + v_surf_component;

            % convert back u,v tp total current speed and direction
            [Uz,Dz] = uv2spddir(u_tot, v_tot, 'to');
            U_total = m_structure(name, xyz, ttt, ['CS z=' num2str(z(zi))], [datenum(tide.Time) Uz], 'Current speed', 1, bins_CS);
            D_total = m_structure(name, xyz, ttt, ['CD z=' num2str(z(zi))], [datenum(tide.Time) Dz], 'Current direction', 1, bins_D);
            U_total.xyz_str = xyz_str;
            D_total.xyz_str = xyz_str;

            % save variables for totals
            save([sdirn '\U_total.mat'],'U_total');
            save([sdirn '\D_total.mat'],'D_total');

        end

    end

end
