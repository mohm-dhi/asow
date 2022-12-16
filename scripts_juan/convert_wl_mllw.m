addpath(genpath('\\USDEN1-STOR.DHI.DK\Projects\41806529\_metocean_scripts\potlab_v2\src'));
addpath(genpath('\\USDEN1-STOR.DHI.DK\Projects\41806529\_metocean_scripts\potlab_v2\res'));

fdir = 'C:\Users\jngz\OneDrive - DHI\2022\41806529-AtlanticShores\Scripts\FileIO\wl\';

% mllw
% mllw = [-0.59,-0.57,-0.57,-0.57,-0.58,-0.59,-0.59]; % dhi
mllw = [-0.55,-0.58,-0.60,-0.61,-0.61,-0.60,-0.63]; %vdatum

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

bins_wl = -2:0.1:2;
ttt = [datenum([1979 01 01 00 00 00]) datenum([2021 12 31 23 00 00]) 30];

for i=1:length(mllw)

    name = ['ASOW' num2str(i)];
    depth = site_depth(i);
    xyz = xyh(i,:);
    xyz_strs = {['(' num2str(xyh(i,1)*-1) 'W;' num2str(xyh(i,2)) 'N;' num2str(xyh(i,3)*-1) 'mMLLW)']};
    xyz_str = xyz_strs{1};

    mkdir([fdir 'MLLW\' name]);

    tidestr = load([fdir name '\' name '_wl.mat']);

    WL = m_structure(name, xyz, ttt, 'WLTotal MLLW', [tidestr.WL.time tidestr.WL.data+abs(mllw(i))], 'WL', 1, bins_wl);
    WLTide = m_structure(name, xyz, ttt, 'WLTide MLLW', [tidestr.WLTde.time tidestr.WLTde.data+abs(mllw(i))], 'WL', 1, bins_wl);
    WLRes = m_structure(name, xyz, ttt, 'WLResidual MLLW', [tidestr.WLRes.time tidestr.WLRes.data], 'WL', 1, bins_wl);

    % edit attributes
    WL.xyz_str = xyz_str;
    WL.vref = 'MLLW';
    WL.unit = 'mMLLW';
    WLTide.xyz_str = xyz_str;
    WLRes.xyz_str = xyz_str;


    % save variables
    save([fdir 'MLLW\' name '\' name '_wl_mllw.mat'],'WL');
    save([fdir 'MLLW\' name '\' name '_wltide_mllw.mat'],'WLTide');
    save([fdir 'MLLW\' name '\' name '_wlres_mllw.mat'],'WLRes');

end