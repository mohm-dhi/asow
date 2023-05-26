% spectra data

model_spec= '\\USDEN1-STOR.DHI.DK\Projects\41806529\02_RAW_MOOD_data\2D Spectra data\Spec0.1_1979-2021_-73.9W39.3N_updated.mat';
obs = 'C:\Workplace\41806529 - Atlantic Shores\Python-Matlab\Spectra_SV2_seg4.mat';
load(model_spec)
load(obs)
obs_datetime = datetime(Time,'ConvertFrom','datenum');

ttt_obs    	= [obs_datetime(1) obs_datetime(end) 10];
ttt_model    	= [obs_datetime(1) obs_datetime(end) 60];
name = 'Spectra';
xyz = [-73.9, 39.3, 0];
legend_model = 'MIKE';
legend_obs = 'Observations';

model_dir = m_structure(name , xyz , ttt_model , legend_model , fname_model ,'CD' , columns_model(2), bins_D);
obs_spd = m_structure(name , xyz , ttt_obs_hd , legend_obs_hd , fname_obs ,'CS' , columns_obs(1), bins_CS);