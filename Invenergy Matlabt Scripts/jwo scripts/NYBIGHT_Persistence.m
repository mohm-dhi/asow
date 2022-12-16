% read in and save data (one time only)
% % E05 tides
% fname='HD E05\US_East_Coast_Water_Level_and_Current_2D_MIKE_21_Hydrodynamic_Model_HD_DHI_-72.71669_39.96928.csv';
% E05_WL = m_structure('E05' , [-72.72 39.97 10] , [datenum(1979,01,01,01,0,0) datenum(2020,12,31,23,0,0) 30] , 'MOOD Water Levels' , fname ,'WL' ,1,[-1.2:0.2:1.2]);
% % E06 tides
% fname='HD E06\US_East_Coast_Water_Level_and_Current_2D_MIKE_21_Hydrodynamic_Model_HD_DHI_-73.42889_39.54677.csv';
% E06_WL = m_structure('E06' , [-73.43 39.55 10] , [datenum([1979 01 01 01 00 00]) datenum([2020 12 31 23 00 00]) 30] , 'MOOD Water Levels' , fname ,'WL' , 1, [-1.2:0.2:1.2]);
% % E05 Currents
% fname='HD E05\US_East_Coast_Water_Level_and_Current_2D_MIKE_21_Hydrodynamic_Model_HD_DHI_-72.71669_39.96928.csv';
% E05_current_speed = m_structure('E05' , [-72.72 39.97 10] , [datenum([1979 01 01 01 00 00]) datenum([2020 12 31 23 00 00]) 30] , 'E05 Total' , fname ,'CS' , 2, [0:0.1:0.6]);
% % E06 Currents
% fname='HD E06\US_East_Coast_Water_Level_and_Current_2D_MIKE_21_Hydrodynamic_Model_HD_DHI_-73.42889_39.54677.csv';
% E06_current_speed = m_structure('E06' , [-72.72 39.97 10] , [datenum([1979 01 01 01 00 00]) datenum([2020 12 31 23 00 00]) 30] , 'E06 Total' , fname ,'CS' , 2, [0:0.1:0.6]);
% % E05 10m winds
% fname='CFSR E05\Global_Met._Parameters_incl._10m_wind_at_0.2_deg._Climate_Forecast_System_Reanalysis_CFSR_NCEP_NOAA_287.2833_39.96928.csv';
% E05_CFSR_10m_speed = m_structure('E05' , [-72.72 39.97 10] , [datenum([1979 01 01 01 00 00]) datenum([2020 12 31 23 00 00]) 60] , 'CFSR 10m' , fname ,'U10' , 1, [0:30]);
% % E06 CFSR 10m winds
% fname='CFSR E06\Global_Met._Parameters_incl._10m_wind_at_0.2_deg._Climate_Forecast_System_Reanalysis_CFSR_NCEP_NOAA_286.5711_39.54677.csv';
% E06_CFSR_10m_speed = m_structure('E06' , [-73.43 39.55 10] , [datenum([1979 01 01 01 00 00]) datenum([2020 12 31 23 00 00]) 60] , 'CFSR 10m' , fname ,'U10' , 1, [0:30]);
% % E05 CFSR 165m
% fname='CFSR Hub Height\E05 CFSR_at_hubHeight_165m.csv';
% E05_CFSR_165m_speed = m_structure('E05' , [-72.72 39.97 10] , [datenum([1979 01 01 01 00 00]) datenum([2020 12 31 23 00 00]) 60] , 'CFSR 165m' , fname ,'U165' , 1, [0:30]);
% % E06 CFSR 165m
% fname='CFSR Hub Height\E06 CFSR_at_hubHeight_165m.csv';
% E06_CFSR_165m_speed = m_structure('E06' , [-73.43 39.55 10] , [datenum([1979 01 01 01 00 00]) datenum([2020 12 31 23 00 00]) 60] , 'CFSR 165m' , fname ,'U165' , 1, [0:30]);
% % E05 Waves
% fname='Waves E05\US_East_Coast_Wave_Parameters_Integrated_MIKE_21_Spectral_Wave_Model_DHI_-72.71669_39.96928.csv';
% E05_waves_Hm0 = m_structure('E05' , [-72.72 39.97 10] , [datenum([1979 01 01 01 00 00]) datenum([2020 12 31 23 00 00]) 60] , 'E05 Total' , fname ,'Hm0' , 1, [0:0.5:12]);
% % E06 Waves
% fname='Waves E06\US_East_Coast_Wave_Parameters_Integrated_MIKE_21_Spectral_Wave_Model_DHI_-73.42889_39.54677.csv';
% E06_waves_Hm0 = m_structure('E06' , [-72.72 39.97 10] , [datenum([1979 01 01 01 00 00]) datenum([2020 12 31 23 00 00]) 60] , 'E06 Total' , fname ,'Hm0' , 1, [0:0.5:12]);
% clear fname
% save PersistenceData
clear all
clc
cd PERSISTENCE2
load PersistenceData
% setup thehresholds and durations
% check call to m_persistence to create EXCEL

% Durations
HDUR = [3 6 9 12 18 24 30 36 48 60 72]; % h

% thresholds
THR_WIND  = [4 6 8 10 12 14]; % m/s
THR_WAVE_Hs  = [0.5 0.75 1.0 1.25 1.5 1.75 2.0 2.25 2.50 2.75 3.0 3.25 3.50 3.75 4.0]; % m
THR_WAVE_Tp  = [4 5 6 7 8]; % m
for PERC=[10 50 90]
    %PERC=50;
    %PERC=90;

    m_persistence(E05_waves_Hm0,'DUR',HDUR,'THR',THR_WAVE_Hs,'P',PERC,'plot_flag',[1]);
    m_persistence(E06_waves_Hm0,'DUR',HDUR,'THR',THR_WAVE_Hs,'P',PERC,'plot_flag',[1]);

    m_persistence(E05_CFSR_10m_speed,'DUR',HDUR,'THR',THR_WIND,'P',PERC,'plot_flag',[1]);
    m_persistence(E06_CFSR_10m_speed,'DUR',HDUR,'THR',THR_WIND,'P',PERC,'plot_flag',[1]);

    m_persistence(E05_waves_Tp,'DUR',HDUR,'THR',THR_WAVE_Tp,'P',PERC,'plot_flag',[1]);
    m_persistence(E06_waves_Tp,'DUR',HDUR,'THR',THR_WAVE_Tp,'P',PERC,'plot_flag',[1]);

end

cd COMBINED
% combined HS and U10 - different limits and only 50% and 90%
THR_WIND  = [6 8 10]; % m/s
THR_WAVE_Hs  = [0.5 1.0 1.5]; % m
for PERC=[50 90]
    for i = 1:length(THR_WIND)
        m_persistence(E05_waves_Hm0,'DUR',HDUR,'THR',THR_WAVE_Hs,'P',PERC,'plot_flag',[1],'Constraint',E05_CFSR_10m_speed,THR_WIND(i),'Below');
        m_persistence(E06_waves_Hm0,'DUR',HDUR,'THR',THR_WAVE_Hs,'P',PERC,'plot_flag',[1],'Constraint',E06_CFSR_10m_speed,THR_WIND(i),'Below');
    end
end
