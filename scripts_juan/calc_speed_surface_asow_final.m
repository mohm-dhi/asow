%%
addpath(genpath('C:\Users\jngz\MATLAB\Projects\potlab\src'));
addpath(genpath('C:\Users\jngz\MATLAB\Projects\potlab\res'));
addpath('C:\DHI\Projects\AtlanticShores\Scripts\asow\scripts_juan\normal_analysis\tidal_analysis');
addpath('C:\DHI\Projects\AtlanticShores\Scripts\asow\scripts_juan\normal_analysis');
addpath('C:\DHI\Projects\AtlanticShores\Scripts\asow\scripts_juan');

%%

fdir = 'C:\DHI\Projects\AtlanticShores\Data\TimeSeries\Structs\';
wdir = 'C:\DHI\Projects\AtlanticShores\Data\TimeSeries\TidalAnalysis\uv\';
no_locs = 7;

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
bins_CS         = 0:0.1:1.5; % current speed
bins_U = -30:0.1:30;
bins_V = -30:0.1:30;
ttt = [datenum([1979 01 15 00 00 00]) datenum([2021 12 31 23 00 00]) 60];

% structs
result = struct;
result_residual = struct;

%%

% for each location
for i=1:no_locs

    cd([wdir 'ASOW' num2str(i)]);

    % load structure
    fname = [fdir 'ASOW' num2str(i) '_all_structs.mat'];
    load(fname);

    % read cfsr wind speed/dir
    cfsr_ws = asow_params.WindSpd_10m;
    cfsr_wd = asow_params.WindDir_10m;

    % read total curr mike speed/dir
    cs_total = asow_params.CurSpd_DepAvg_Total;
    cd_total = asow_params.CurDir_DepAvg_Total;

    % read tide curr mike speed/dir
    cs_tide = asow_params.CurSpd_DepAvg_Tide;
    cd_tide = asow_params.CurDir_DepAvg_Tide;

    % read atmospheric curr mike speed/dir
    cs_res = asow_params.CurSpd_DepAvg_Res;
    cd_res = asow_params.CurDir_DepAvg_Res;

    % descr
    name = ['ASOW' num2str(i)];
    result.(name) = [];
    result_residual.(name) = [];
    depth = site_depth(i);
    xyz = xyh(i,:);
    xyz_strs = {['(' num2str(xyh(i,1)*-1) 'W;' num2str(xyh(i,2)*-1) 'N;' num2str(xyh(i,3)*-1) 'mMSL)']};
    xyz_str = xyz_strs{1};

    % calculate offsets for desired depths
    % fr = ratio of depth (ie 0.75 = 75% upwards from depth; 1 = surface)
    fr = [1,0.75,0.5,0.25,0.05];

    % z = distance from top, 0 is surface
    z = (depth.*fr) - depth;

    % calculate u_wind
    u_wind = 0.01 * cfsr_ws.data;

    % apply depth correction at all levels
    u_wind_z = u_wind_correction(u_wind, z, 20);

    % integrate for depth-averaged speed of wind
    zb = depth + z; % integration depths
    total_depth = site_depth(i); % total depth at site
    u_wind_bar = ((1/total_depth)*trapz(fliplr(zb),flipud(u_wind_z),1))'; % integrated speed
    u_wind_dir = mod(cfsr_wd.data+180,360); % direction same as cfsr, towards

    % calculate u_atm
    cs_atm = cs_res.data - u_wind_bar;
    cd_atm = cd_res.data;

    % depth levels
    % level 1 for surface
    % level 5 for bottom
    for zi = 1:length(z)

        % make directory for level
        sdirn = ['level_' num2str(zi)];
        mkdir(sdirn);

        % zb = 0 at bottom
        % zb = depth at surface
        zb = depth + z(zi);

        disp(["zb" num2str(zb) "z" num2str(z(zi))])

        if zb>0

            % calculate ubar_tide at zb
            if depth/2<zb & zb<depth
                ubar_tide_zb = cs_tide.data* 1.07;
            else
                ubar_tide_zb = cs_tide.data*(zb/(0.32*depth))^(1/7);
            end

            % calculate ubar_atm at zb
            if depth/2<zb & zb<depth
                ubar_atm_zb = cs_atm* 1.07;
            else
                ubar_atm_zb = cs_atm*(zb/(0.32*depth))^(1/7);
            end

            % calculate u residual (ubar_atm_zb + u_wind_zb)
            [u1,v1] = spddir2uv(ubar_atm_zb,cd_atm);
            [u2,v2] = spddir2uv(u_wind_z(zi,:)',u_wind_dir);

            u_res_z_uv = [u1+u2,v1+v2];
            [u_res_z_sp,u_res_z_dir] = uv2spddir(u_res_z_uv(:,1),u_res_z_uv(:,2),'to');

            % calculate u total (tide + residual)
            [u3,v3] = spddir2uv(ubar_tide_zb,cd_tide.data);

            u_total_z_uv = [u_res_z_uv(:,1)+u3,u_res_z_uv(:,2)+v3];
            [u_total_z_sp,u_total_z_dir] = uv2spddir(u_total_z_uv(:,1),u_total_z_uv(:,2),'to');

            % make m_structure for tidal current at depth
            U_tide = m_structure(name, xyz, ttt,'CS_{Tide}',[datenum(cs_tide.time) ubar_tide_zb cd_tide.data], 'Current_speed_tide', 1, bins_CS);
            D_tide = m_structure(name, xyz, ttt,'CD_{Tide}',[datenum(cs_tide.time) ubar_tide_zb cd_tide.data], 'Current_direction_tide', 2, bins_D);

            save([sdirn '\U_tide.mat'],'U_tide');
            save([sdirn '\D_tide.mat'],'D_tide');

            % make structures for residual current
            U_res = m_structure(name, xyz, ttt, ['CS_{Residual} z=' num2str(z(zi))], [datenum(cs_tide.time) u_res_z_sp], 'Current_speed_residual', 1, bins_CS);
            D_res = m_structure(name, xyz, ttt, ['CD_{Residual} z=' num2str(z(zi))], [datenum(cs_tide.time) u_res_z_dir], 'Current_direction_residual', 1, bins_D);
            U_res.xyz_str = xyz_str;
            D_res.xyz_str = xyz_str;

            % save variables for residuals
            save([sdirn '\U_res.mat'],'U_res');
            save([sdirn '\D_res.mat'],'D_res');

            U_total = m_structure(name, xyz, ttt, ['CS z=' num2str(z(zi))], [datenum(cs_tide.time) u_total_z_sp], 'Current speed', 1, bins_CS);
            D_total = m_structure(name, xyz, ttt, ['CD z=' num2str(z(zi))], [datenum(cs_tide.time) u_total_z_dir], 'Current direction', 1, bins_D);
            U_total.xyz_str = xyz_str;
            D_total.xyz_str = xyz_str;

            % save variables for totals
            save([sdirn '\U_total.mat'],'U_total');
            save([sdirn '\D_total.mat'],'D_total');

        end

    end

end
