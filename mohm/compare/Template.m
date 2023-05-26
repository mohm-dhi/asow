load cmp_data.mat

if contains(type, 'current')
    model_spd = m_structure(name , xyz , ttt_model_hd , legend_model_hd , fname_model ,'CS' , columns_model(1), bins_CS);
    model_dir = m_structure(name , xyz , ttt_model_hd , legend_model_hd , fname_model ,'CD' , columns_model(2), bins_D);
    obs_spd = m_structure(name , xyz , ttt_obs , legend_obs_hd , fname_obs ,'CS' , columns_obs(1), bins_CS);
    obs_dir = m_structure(name , xyz , ttt_obs , legend_obs_hd , fname_obs ,'CD' , columns_obs(2), bins_D);
    
    model_spd.ascii = '.xlsx';
    model_dir.ascii = '.xlsx';
    obs_spd.ascii = '.xlsx';
    obs_dir.ascii = '.xlsx';
    obs_spd.dfs_info = ' ';
    obs_dir.dfs_info = ' ';
    % Timeseries
    
    m_timeseries([obs_spd model_spd ],[ obs_dir model_dir]);
    m_compare(obs_spd,model_spd, 'dir',obs_dir, model_dir, 'bin', [0 0.05:0.05:0.4 0.9]);
    m_compare(obs_dir, model_dir, 'dir',obs_dir, model_dir);
elseif contains(type, 'wave')
    model_Hm0 = m_structure(name , xyz , ttt_model_sw , legend_model_sw , fname_model ,'Hm0' , columns_model(1), bins_H);
    model_Tp = m_structure(name , xyz , ttt_model_sw , legend_model_sw , fname_model ,'Tp' , columns_model(2), bins_T);
    model_MWD = m_structure(name , xyz , ttt_model_sw , legend_model_sw , fname_model ,'MWD' , columns_model(3), bins_D);
    obs_Hm0 = m_structure(name , xyz , ttt_obs , legend_obs_sw , fname_obs ,'Hm0' , columns_obs(1), bins_H);
    obs_Tp = m_structure(name , xyz , ttt_obs , legend_obs_sw , fname_obs ,'Tp' , columns_obs(2), bins_T);
    obs_MWD = m_structure(name , xyz , ttt_obs , legend_obs_sw , fname_obs ,'MWD' , columns_obs(3), bins_D);
    
	obs_Hm0.data = obs_Hm0.data.^2;
	obs_Hm0=m_resample(obs_Hm0,'dt_int',60,'dt_avg',180);
	obs_Hm0.data = sqrt(obs_Hm0.data);
    							   
    model_Hm0.ascii = '.xlsx';
    model_Tp.ascii = '.xlsx';
    model_MWD.ascii = '.xlsx';
    obs_Hm0.ascii = '.xlsx';
    obs_Tp.ascii = '.xlsx';
    obs_MWD.ascii = '.xlsx';
    obs_Hm0.dfs_info = ' ';
    obs_Tp.dfs_info = ' ';
    obs_MWD.dfs_info = ' ';
    % Timeseries
    
    m_timeseries([ obs_Hm0 model_Hm0 ],[ obs_Tp model_Tp], [ obs_MWD model_MWD]);
    m_compare(obs_Hm0, model_Hm0,'dir',obs_MWD, model_MWD,  'bin', [0 0:0.5:6 7])
    m_compare(obs_Tp, model_Tp, 'dir',obs_MWD, model_MWD)
    m_compare(obs_MWD, model_MWD, 'dir',obs_MWD, model_MWD)
elseif contains (type, 'wl')
    model_wl = m_structure(name , xyz , ttt_model , legend_model_hd , fname_model ,'WL' , columns_model(1), bins_WL);
    obs_wl = m_structure(name , xyz , ttt_obs , legend_obs_hd , fname_obs ,'WL' , columns_obs(1), bins_WL);
    model_wl.data = model_wl.data - nanmean(model_wl.data);
    obs_wl.data = obs_wl.data - nanmean(obs_wl.data);
    model_wl.ascii = '.xlsx';
    obs_wl.ascii = '.xlsx';
    obs_wl.dfs_info = ' ';
    
    % Timeseries
    m_timeseries([obs_wl model_wl ]);
    m_compare(obs_wl, model_wl)
    
    
elseif contains(type, 'wind')
    name_saved = name;
    name = [name ' at 20m'];
    model_spd = m_structure(name , xyz , ttt_model , legend_model_wind , fname_model ,'WS' , columns_model(1), bins_WS);
    U20 = U102Uz(model_spd.data,20,'power',0.111);
    model_spd.data = U20;
    model_dir = m_structure(name , xyz , ttt_model , legend_model_wind , fname_model ,'WD' , columns_model(2), bins_D);
    obs_spd = m_structure(name , xyz , ttt_obs , legend_obs_wind , fname_obs ,'WS' , columns_obs(1), bins_WS);
    obs_dir = m_structure(name , xyz , ttt_obs , legend_obs_wind , fname_obs ,'WD' , columns_obs(2), bins_D);
    
    [obs_u, obs_v] = spddir2uv(obs_spd, obs_dir, 'invert');    
    obs_u = m_resample(obs_u,'dt_avg',120);
    obs_v = m_resample(obs_v,'dt_avg',120);
    [obs_spd, obs_dir] = uv2spddir(obs_u, obs_v, 'from');															   
    model_spd.ascii = '.xlsx';
    model_dir.ascii = '.xlsx';
    obs_spd.ascii = '.xlsx';
    obs_dir.ascii = '.xlsx';
    obs_spd.dfs_info = ' ';
    obs_dir.dfs_info = ' ';
    % Timeseries
    
    m_timeseries([obs_spd model_spd ],[obs_dir model_dir ]);
    m_compare(obs_spd, model_spd,'dir',obs_dir ,model_dir, 'bin', [0:3:33]);
    
    
    name = [name_saved ' at 40m'];
    model_spd = m_structure(name , xyz , ttt_model , legend_model_wind , fname_model ,'WS' , columns_model(1), bins_WS);
    U20 = U102Uz(model_spd.data,40,'power',0.111);
    model_spd.data = U20;
    model_dir = m_structure(name , xyz , ttt_model , legend_model_wind , fname_model ,'WD' , columns_model(2), bins_D);
    obs_spd = m_structure(name , xyz , ttt_obs , legend_obs_wind , fname_obs ,'WS' , columns_obs(3), bins_WS);
    obs_dir = m_structure(name , xyz , ttt_obs , legend_obs_wind , fname_obs ,'WD' , columns_obs(4), bins_D);
    
    [obs_u, obs_v] = spddir2uv(obs_spd, obs_dir, 'invert');    
    obs_u = m_resample(obs_u,'dt_avg',120);
    obs_v = m_resample(obs_v,'dt_avg',120);
    [obs_spd, obs_dir] = uv2spddir(obs_u, obs_v, 'from');															   
    model_spd.ascii = '.xlsx';
    model_dir.ascii = '.xlsx';
    obs_spd.ascii = '.xlsx';
    obs_dir.ascii = '.xlsx';
    obs_spd.dfs_info = ' ';
    obs_dir.dfs_info = ' ';
    % Timeseries
    
    m_timeseries([obs_spd model_spd ],[obs_dir model_dir ]);
    m_compare(obs_spd, model_spd,'dir',obs_dir ,model_dir, 'bin', [0:3:33]);
    
    

    name = [name_saved ' at 60m'];
    model_spd = m_structure(name , xyz , ttt_model , legend_model_wind , fname_model ,'WS' , columns_model(1), bins_WS);
    U20 = U102Uz(model_spd.data,60,'power',0.111);
    model_spd.data = U20;
    model_dir = m_structure(name , xyz , ttt_model , legend_model_wind , fname_model ,'WD' , columns_model(2), bins_D);
    obs_spd = m_structure(name , xyz , ttt_obs , legend_obs_wind , fname_obs ,'WS' , columns_obs(5), bins_WS);
    obs_dir = m_structure(name , xyz , ttt_obs , legend_obs_wind , fname_obs ,'WD' , columns_obs(6), bins_D);
    
    [obs_u, obs_v] = spddir2uv(obs_spd, obs_dir, 'invert');    
    obs_u = m_resample(obs_u,'dt_avg',120);
    obs_v = m_resample(obs_v,'dt_avg',120);
    [obs_spd, obs_dir] = uv2spddir(obs_u, obs_v, 'from');															   
    model_spd.ascii = '.xlsx';
    model_dir.ascii = '.xlsx';
    obs_spd.ascii = '.xlsx';
    obs_dir.ascii = '.xlsx';
    obs_spd.dfs_info = ' ';
    obs_dir.dfs_info = ' ';
    % Timeseries
    
    m_timeseries([obs_spd model_spd ],[obs_dir model_dir ]);
    m_compare(obs_spd, model_spd,'dir',obs_dir ,model_dir, 'bin', [0:3:33]);
    
    
    name = [name_saved ' at 80m'];
    model_spd = m_structure(name , xyz , ttt_model , legend_model_wind , fname_model ,'WS' , columns_model(1), bins_WS);
    U20 = U102Uz(model_spd.data,80,'power',0.111);
    model_spd.data = U20;
    model_dir = m_structure(name , xyz , ttt_model , legend_model_wind , fname_model ,'WD' , columns_model(2), bins_D);
    obs_spd = m_structure(name , xyz , ttt_obs , legend_obs_wind , fname_obs ,'WS' , columns_obs(7), bins_WS);
    obs_dir = m_structure(name , xyz , ttt_obs , legend_obs_wind , fname_obs ,'WD' , columns_obs(8), bins_D);
    [obs_u, obs_v] = spddir2uv(obs_spd, obs_dir, 'invert');    
    obs_u = m_resample(obs_u,'dt_avg',120);
    obs_v = m_resample(obs_v,'dt_avg',120);
    [obs_spd, obs_dir] = uv2spddir(obs_u, obs_v, 'from');															   
    
    model_spd.ascii = '.xlsx';
    model_dir.ascii = '.xlsx';
    obs_spd.ascii = '.xlsx';
    obs_dir.ascii = '.xlsx';
    obs_spd.dfs_info = ' ';
    obs_dir.dfs_info = ' ';
    % Timeseries
    
    m_timeseries([obs_spd model_spd ],[obs_dir model_dir ]);
    m_compare(obs_spd, model_spd,'dir',obs_dir ,model_dir, 'bin', [0:3:33]);
    m_compare(obs_dir ,model_dir,'dir',obs_dir ,model_dir)
    
end

