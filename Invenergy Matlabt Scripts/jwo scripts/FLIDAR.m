    clc
    clear all
    % read in just the 160m FLIDAR data WS, WS std and TI at 160m for TI_model call - no NANs
    fname='NYSERDAFLiDAR\Data\E05_Hudson_North_10_min_avg_20190812_20210919_160m.csv';
    E05_WS_160m_FLIDAR = m_structure('E05' , [-72.72 39.97 160] , [datenum(2019,08,12,0,0,0) datenum(2021,9,19,23,50,0) 10] , 'FLIDAR Speed' , fname ,'U100' ,1,[0:2:30]);
    E05_WD_160m_FLIDAR = m_structure('E05' , [-72.72 39.97 160] , [datenum(2019,08,12,0,0,0) datenum(2021,9,19,23,50,0) 10] , 'FLIDAR Direction' , fname ,'D100' ,3,[0:30:360]);

    E05_WS_Std_160m_FLIDAR = m_structure('E05' , [-72.72 39.97 160] , [datenum(2019,08,12,0,0,0) datenum(2021,9,19,23,50,0) 10] , 'FLIDAR Std' , fname ,'WSSTD160' ,2,[0:.1:2]);
    E05_WS_TI_160m_FLIDAR = m_structure('E05' , [-72.72 39.97 160] , [datenum(2019,08,12,0,0,0) datenum(2021,9,19,23,50,0) 10] , 'FLIDAR TI' , fname ,'TI160',4,[0:.1:2]);
    % replace TI in raw data with std/V ??
    E05_WS_TI_160m_FLIDAR.data=E05_WS_Std_160m_FLIDAR.data./E05_WS_160m_FLIDAR.data;

    fname='..\NYSERDAFLiDAR\Data\E06_Hudson_South_10_min_avg_20190904_20220123_160m.csv';
    E06_WS_160m_FLIDAR = m_structure('E06' , [-72.72 39.97 160] , [datenum(2019,9,4,0,0,0) datenum(2022,1,23,23,50,0) 10] , 'FLIDAR WS' , fname ,'U100' ,1,[0:2:30]);
    E06_WS_Std_160m_FLIDAR = m_structure('E06' , [-72.72 39.97 160] , [datenum(2019,9,4,0,0,0) datenum(2022,1,23,23,50,0) 10] , 'FLIDAR Std' , fname ,'WSSTD160' ,2,[0:.1:2]);
    E06_WD_160m_FLIDAR = m_structure('E06' , [-72.72 39.97 160] , [datenum(2019,9,4,0,0,0) datenum(2022,1,23,23,50,0) 10] , 'FLIDAR Direction' , fname ,'D100' ,3,[0:30:360]);
    E06_WS_TI_160m_FLIDAR = m_structure('E06' , [-72.72 39.97 160] , [datenum(2019,9,4,0,0,0) datenum(2022,1,23,23,50,0) 10] , 'FLIDAR TI' , fname ,'TI160' ,4,[0:.1:2]);
    % replace TI in raw data with std/V ??
    E06_WS_TI_160m_FLIDAR.data=E06_WS_Std_160m_FLIDAR.data./E06_WS_160m_FLIDAR.data;


E05_WD_160m_FLIDAR.legend='FLIDAR at 160m';
E05_WD_160m_FLIDAR.label='D_{100}';
% do scatter plot
m_scatter(E05_WD_160m_FLIDAR,E05_WS_160m_FLIDAR,'ci',0:0.05:.25)
