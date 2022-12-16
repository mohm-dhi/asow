load cmp_data.mat

SW_WS = m_structure(name , xyz , ttt_model , legend_model , fname_model ,'WS' , columns_model(1), bins_WS);
SW_WS.T = [1 5 10 50 100 500 1000 10000];
sensitivity = 1;
SW_WS.fontsize = 6;

if sensitivity
    m_timeseries(SW_WS)
    [Ip_ALL,ID,Nyr] = m_peak(SW_WS, 'AMP',2, 36, 0.7, 'plotflag',1);
    m_extreme_sensitivity(SW_WS)
else
    opt_CS.EVdist           = 3;
    opt_CS.EVtype           = 'AAP';
    opt_CS.EVcrit           = 3;
    opt_CS.N_BootStrap      = 100;
    opt_CS.estimationmethod = 'LS';
    opt_CS.plotflag         = [1]; %6 [1,2,4]
    opt_CS.EVadju           = 0;
    opt_CS.ndec             = 1;
    opt_CS.constfac         = 1e2;
    Opt_Hm0.optstr           = 'directional';
%     m_extreme(m_subseries(SW_CS,'directional',SW_CS),opt_CS);% directional
%     m_extreme(m_subseries(SW_CS,'monthly'),opt_CS); %directional
%     Opt_CS.optstr           = 'omni';
    m_extreme(SW_CS ,opt_CS); %omnidirectional
end