
clear all
close all
warning off
% MIKE stuff
NET.addAssembly('DHI.Mike.Install');
import DHI.Mike.Install.*;
DHI.Mike.Install.MikeImport.SetupLatest({DHI.Mike.Install.MikeProducts.MikeCore});
% load parameters from readNYBIGHT script
load NYBIGHT

% plots for NYBIGHT report

%% ice accretion
cd ICE
% ice risk criteria - use TI call as I've edited it to leave figure open
% add ice risk lines onto scatter plot and save to jpeg
E05_CFSR_10m_SST.legend='';
E05_CFSR_10m_temperature.legend='2 m';
E05_CFSR_10m_SST.bins=[-15:3:30];
E05_CFSR_10m_temperature.bins=[-15:3:30];
m_scatter_TI(E05_CFSR_10m_temperature,E05_CFSR_10m_SST,'ci',-15:3:30,'density')
fprintf(1,'%s %g\n','%Time air temp. less than -2^oC (E05)',100*length(find(E05_CFSR_10m_temperature.data<-2))/length(E05_CFSR_10m_temperature.data));
fprintf(1,'%s %g\n','%Time air temp. less than -2^oC (E06)',100*length(find(E06_CFSR_10m_temperature.data<-2))/length(E06_CFSR_10m_temperature.data));
fprintf(1,'%s %g\n','%Time SST less than 8^oC (E05)',100*length(find(E05_CFSR_10m_SST.data<8))/length(E05_CFSR_10m_SST.data));
fprintf(1,'%s %g\n','%Time SST less than 8^oC (E06)',100*length(find(E06_CFSR_10m_SST.data<8))/length(E06_CFSR_10m_SST.data));

IsRisky=100*(length(find(E05_CFSR_10m_temperature.data<-2 & E05_CFSR_10m_SST.data<8))/length(E05_CFSR_10m_temperature.data));
hold on
plot([-2 -2 -15],[-15 8 8],'k--','lineWidth',1,'DisplayName','Ice Risk Window');
text(-10,-4,'Risk of ice');
text(-10,-6,[num2str(IsRisky,2), '%']);
text(3,-10,'No risk of ice');
print -djpeg E05Ice_risk.jpg;close

clear IsRisky
m_scatter_TI(E06_CFSR_10m_temperature,E06_CFSR_10m_SST,'ci',-15:3:30,'density')
E06_CFSR_10m_SST.legend='';
E06_CFSR_10m_temperature.legend='2 m';
E06_CFSR_10m_SST.bins=[-15:3:30];
E06_CFSR_10m_temperature.bins=[-15:3:30];
IsRisky=100*(length(find(E06_CFSR_10m_temperature.data<-2 & E06_CFSR_10m_SST.data<8))/length(E06_CFSR_10m_temperature.data));
hold on
plot([-2 -2 -15],[-15 8 8],'k--','lineWidth',1,'DisplayName','Ice Risk Window');
text(-10,-4,'Risk of ice');
text(-10,-6,[num2str(IsRisky,2), '%']);
text(3,-10,'No risk of ice');
print -djpeg E06Ice_risk.jpg;close

fprintf(1,'%s %g %g\n','E05 min/max SST',min(E05_CFSR_10m_SST.data),max(E05_CFSR_10m_SST.data))
fprintf(1,'%s %g %g\n','E06 min/max SST',min(E06_CFSR_10m_SST.data),max(E06_CFSR_10m_SST.data))
fprintf(1,'%s %g %g\n','E05 min/max air temperature',min(E05_CFSR_10m_temperature.data),max(E05_CFSR_10m_temperature.data))
fprintf(1,'%s %g %g\n','E06 min/max air temperature',min(E06_CFSR_10m_temperature.data),max(E06_CFSR_10m_temperature.data))


cd ..

%% plot 10m wind roses
cd ROSES
m_scatter(E05_CFSR_10m_direction,E05_CFSR_10m_speed)
%m_rose_plot(E05_CFSR_10m_direction.data,E05_CFSR_10m_speed.data,'Ay',[0:20:360],'di',[2:18],'ci',[5:5:15])
%print -djpeg E05_10m_wind_rose.jpg;close
m_scatter(E06_CFSR_10m_direction,E06_CFSR_10m_speed)
%m_rose_plot(E06_CFSR_10m_direction.data,E06_CFSR_10m_speed.data,'Ay',[0:20:360],'di',[2:18],'ci',[5:5:15])
%print -djpeg E06_10m_wind_rose.jpg;close
cd ..
%% plot 100 wind roses - not needed for report
% m_rose_plot(E05_CFSR_100m_direction.data,E05_CFSR_100m_speed.data,'Ay',[0:20:360],'di',[2:18],'ci',[5:5:15])
% print -djpeg E05_100m_wind_rose.jpg;close
% m_rose_plot(E06_CFSR_100m_direction.data,E06_CFSR_100m_speed.data,'Ay',[0:20:360],'di',[2:18],'ci',[5:5:15])
% print -djpeg E06_100m_wind_rose.jpg;close

%% CFSR 10m + histograms plots
% with Weibull distribution
% figures with 2xWeibull
cd PROBABILITY
m_probability_Weibull(E05_CFSR_10m_speed);
% directional
[outputStat,outputPDF]=m_probability_Weibull(E05_CFSR_10m_speed,'directional',E05_CFSR_10m_direction)
m_probability_Weibull(E06_CFSR_10m_speed);
[outputStat,outputPDF]=m_probability_Weibull(E06_CFSR_10m_speed,'directional',E06_CFSR_10m_direction)


%% CFSR  165m + histograms plots
% change label and units - not right in the m_structure read!
E05_CFSR_165m_speed.legend='CFSR at 165m';
E05_CFSR_165m_speed.label='U_{165}';
E05_CFSR_165m_speed.unit='m/s';
m_probability_Weibull(E05_CFSR_165m_speed);
[outputStat,outputPDF]=m_probability_Weibull(E05_CFSR_165m_speed,'directional',E05_CFSR_165m_direction)
E06_CFSR_165m_speed.legend='CFSR at 165m';
E06_CFSR_165m_speed.label='U_{165}';
E06_CFSR_165m_speed.unit='m/s';
m_probability_Weibull(E06_CFSR_165m_speed);
[outputStat,outputPDF]=m_probability_Weibull(E06_CFSR_165m_speed,'directional',E06_CFSR_165m_direction)
cd ..
%% 165m wind roses
cd ROSES
E05_CFSR_165m_direction.legend='CFSR at 165m';
E05_CFSR_165m_direction.label='D_{165}';
m_scatter(E05_CFSR_165m_direction,E05_CFSR_165m_speed,'di',[2:18],'ci',[5:5:15])
%m_rose_plot(E05_CFSR_165m_direction.data,E05_CFSR_165m_speed.data,'Ay',[0:20:360],'di',[2:18],'ci',[5:5:15])
%print -djpeg E05_165m_wind_rose.jpg;close
E06_CFSR_165m_direction.legend='CFSR at 165m';
E06_CFSR_165m_direction.label='D_{165}';
m_scatter(E06_CFSR_165m_direction,E06_CFSR_165m_speed,'di',[2:18],'ci',[5:5:15])
%m_rose_plot(E06_CFSR_165m_direction.data,E06_CFSR_165m_speed.data,'Ay',[0:20:360],'di',[2:18],'ci',[5:5:15])
%print -djpeg E06_165m_wind_rose.jpg;close
cd ..
% look at profiles - not used for report but of interest
plot([mean(E05_CFSR_10m_speed.data) mean(E05_CFSR_100m_speed.data) mean(E05_CFSR_165m_speed.data) ],[10 100 165]);
hold on 
plot([mean(E06_CFSR_10m_speed.data) mean(E06_CFSR_100m_speed.data) mean(E06_CFSR_165m_speed.data)],[10 100 165]);
grid on
title('CFSR 10m, 100m and extrapolated to 165m)')
legend('E05','E06','Location','NorthWest')
print -djpeg meancfsr_winds.jpg;close

%% create contour plots of wind speed + std by month/sector
% split into month and speed bins
% both 100m and 165m in here but only use 165m in report
cd WINDCONTOUR
tt=datevec(E05_CFSR_165m_speed.time);
for I=1:12
    for K=0:30:330
        thismonthsector=find(tt(:,2)==I & E05_CFSR_165m_direction.data>=K & E05_CFSR_165m_direction.data<K+30);
        E05Speed165m(I,K/30+1)=mean(E05_CFSR_165m_speed.data(thismonthsector));
        E05Std165m(I,K/30+1)=std(E05_CFSR_165m_speed.data(thismonthsector));
        clear thismonthsector
    end
end
contour([1:12],[0:30:330],E05Speed165m,'LevelList',[4:.5:12],'ShowText','on','LineWidth',2)
ylabel('Sector');
set(gca,'xticklabel', {['Jan'], ['Feb'], ['Mar'],['Apr'],['May'], ['Jun'],['Jul'], ['Aug'], ['Sep'], ['Oct'],['Nov'], ['Dec'], ''});
print -djpeg E05_165m_Sector_speed.jpg;close

contour([1:12],[0:30:330],E05Std165m,'LevelList',[2:0.2:5.2],'ShowText','on','LineWidth',2)
ylabel('Sector');
set(gca,'xticklabel', {['Jan'], ['Feb'], ['Mar'],['Apr'],['May'], ['Jun'],['Jul'], ['Aug'], ['Sep'], ['Oct'],['Nov'], ['Dec'], ''});
print -djpeg E05_165m_Sector_std.jpg;close


clear tt
tt=datevec(E06_CFSR_165m_speed.time);
for I=1:12
    for K=0:30:330
        thismonthsector=find(tt(:,2)==I & E06_CFSR_165m_direction.data>=K & E06_CFSR_165m_direction.data<K+30);
        E06Speed165m(I,K/30+1)=mean(E06_CFSR_165m_speed.data(thismonthsector));
        E06Std165m(I,K/30+1)=std(E06_CFSR_165m_speed.data(thismonthsector));
        clear thismonthsector
    end
end
contour([1:12],[0:30:330],E06Speed165m,'LevelList',[4:.5:12],'ShowText','on','LineWidth',2)
ylabel('Sector');
set(gca,'xticklabel', {['Jan'], ['Feb'], ['Mar'],['Apr'],['May'], ['Jun'],['Jul'], ['Aug'], ['Sep'], ['Oct'],['Nov'], ['Dec'], ''});
print -djpeg E06_165m_Sector_speed.jpg;close

contour([1:12],[0:30:330],E06Std165m,'LevelList',[2:0.2:5.2],'ShowText','on','LineWidth',2)
ylabel('Sector');
set(gca,'xticklabel', {['Jan'], ['Feb'], ['Mar'],['Apr'],['May'], ['Jun'],['Jul'], ['Aug'], ['Sep'], ['Oct'],['Nov'], ['Dec'], ''});
print -djpeg E06_165m_Sector_std.jpg;close


tt=datevec(E05_CFSR_100m_speed.time);
for I=1:12
    for K=0:30:330
        thismonthsector=find(tt(:,2)==I & E05_CFSR_100m_direction.data>=K & E05_CFSR_100m_direction.data<K+30);
        E05Speed100m(I,K/30+1)=mean(E05_CFSR_100m_speed.data(thismonthsector));
        E05Std100m(I,K/30+1)=std(E05_CFSR_100m_speed.data(thismonthsector));
        clear thismonthsector
    end
end
contour([1:12],[0:30:330],E05Speed100m,'LevelList',[4:.5:12],'ShowText','on','LineWidth',2)
ylabel('Sector');
set(gca,'xticklabel', {['Jan'], ['Feb'], ['Mar'],['Apr'],['May'], ['Jun'],['Jul'], ['Aug'], ['Sep'], ['Oct'],['Nov'], ['Dec'], ''});
print -djpeg E05_100m_Sector_speed.jpg;close

contour([1:12],[0:30:330],E05Std100m,'LevelList',[2:0.2:5.2],'ShowText','on','LineWidth',2)
ylabel('Sector');
set(gca,'xticklabel', {['Jan'], ['Feb'], ['Mar'],['Apr'],['May'], ['Jun'],['Jul'], ['Aug'], ['Sep'], ['Oct'],['Nov'], ['Dec'], ''});
print -djpeg E05_100m_Sector_std.jpg;close


tt=datevec(E06_CFSR_100m_speed.time);
for I=1:12
    for K=0:30:330
        thismonthsector=find(tt(:,2)==I & E06_CFSR_100m_direction.data>=K & E06_CFSR_100m_direction.data<K+30);
        E06Speed100m(I,K/30+1)=mean(E06_CFSR_100m_speed.data(thismonthsector));
        E06Std100m(I,K/30+1)=std(E06_CFSR_100m_speed.data(thismonthsector));
        clear thismonthsector
    end
end
contour([1:12],[0:30:330],E06Speed100m,'LevelList',[4:.5:12],'ShowText','on','LineWidth',2)
ylabel('Sector');
set(gca,'xticklabel', {['Jan'], ['Feb'], ['Mar'],['Apr'],['May'], ['Jun'],['Jul'], ['Aug'], ['Sep'], ['Oct'],['Nov'], ['Dec'], ''});
print -djpeg E06_100m_Sector_speed.jpg;close

contour([1:12],[0:30:330],E06Std100m,'LevelList',[2:0.2:5.2],'ShowText','on','LineWidth',2)
ylabel('Sector');
set(gca,'xticklabel', {['Jan'], ['Feb'], ['Mar'],['Apr'],['May'], ['Jun'],['Jul'], ['Aug'], ['Sep'], ['Oct'],['Nov'], ['Dec'], ''});
print -djpeg E06_100m_Sector_std.jpg;close
cd ..

%% TI plots
cd TURBULENCE
% check on CNT and H0 TI_model script - used values from m_scatter_turb.m
% load in the apprpraite windspped and standard deviation data
WS=E05_WS_160m_FLIDAR.data';
TI=E05_WS_TI_160m_FLIDAR.data';
m_rose_plot(E05_WD_160m_FLIDAR.data,E05_WS_TI_160m_FLIDAR.data,'Ay',[0:30:360],'ci',0:5:15,'di',[0 0.01 0.02 0.05 0.1 0.2 0.3],'labtitle',['E05 TI at 160m'] )
print -djpeg E05TIrose.jpg;close
% call TI_models
[TI_offshore,TI_90_offshore,TI_IEC_NTM_A,TI_IEC_NTM_B,TI_IEC_NTM_C,TI_DNV_NTM_OA,TI_DNV_NTM_OB,TI_DNV_NTM_OC] = TI_models(165,WS,2.4,24.2);

% split into month and time of day
for I=1:32
    thisspeedbin=find(WS>I-1 & WS<I);
    E05TISPEED(I)=mean(TI(thisspeedbin));
%    E05TISTDSPEED(I)=std(TI(thisspeedbin))/sqrt(length(thisspeedbin));
    E05TISTDSPEED(I)=std(TI(thisspeedbin));
    E05TI90SPEED(I)=prctile(TI(thisspeedbin),90);
    clear  thisspeedbin
end

% do coloued scatted plot
m_scatter_TI(E05_WS_160m_FLIDAR,E05_WS_TI_160m_FLIDAR,'ci',0:0.05:.25,'density','weights',ones(99999,1));
legend('Location','NorthEast')
hold on
% workaround to plot TI estimates
errorbar(E05TISPEED,E05TISTDSPEED,'r','lineWidth',1,'DisplayName','TI')
plot([1:32],E05TI90SPEED,'r-','lineWidth',1,'DisplayName','TI 90')
% sort the estiamtes into order so plot looks like a single line
[B,I]=sort(WS);
plot(WS(I),TI_offshore(I),'lineWidth',1,'DisplayName','TI offhsore');
plot(WS(I),TI_90_offshore(I),'lineWidth',1,'DisplayName','TI offshore 90');
plot(WS(I),TI_DNV_NTM_OA(I),'lineWidth',1,'DisplayName','DNV NTM OA');
plot(WS(I),TI_DNV_NTM_OB(I),'lineWidth',1,'DisplayName','DNV NTM OB');
h1=plot(WS(I),TI_DNV_NTM_OC(I),'lineWidth',1,'DisplayName','DNV NTM OC');
set(gca,'Ylim',[0 .7]);
print -djpeg E05TIscatter.jpg;close

% repeat for E06
clear WS TI TI_offshore TI_90_offshore TI_IEC_NTM_A TI_IEC_NTM_B TI_IEC_NTM_C TI_DNV_NTM_OA TI_DNV_NTM_OB TI_DNV_NTM_OC
WS=E06_WS_160m_FLIDAR.data';
TI=E06_WS_TI_160m_FLIDAR.data';
m_rose_plot(E06_WD_160m_FLIDAR.data,E06_WS_TI_160m_FLIDAR.data,'Ay',[0:30:360],'ci',0:5:15,'di',[0 0.01 0.02 0.05 0.1 0.2 0.3],'labtitle',['E06 TI at 160m'] )
print -djpeg E06TIrose.jpg;close

[TI_offshore,TI_90_offshore,TI_IEC_NTM_A,TI_IEC_NTM_B,TI_IEC_NTM_C,TI_DNV_NTM_OA,TI_DNV_NTM_OB,TI_DNV_NTM_OC] = TI_models(165,WS,2.4,24.2);
% workaround to plot estimates
for I=1:32
    thisspeedbin=find(WS>I-1 & WS<I);
    E06TISPEED(I)=mean(TI(thisspeedbin));
%    E06TISTDSPEED(I)=std(TI(thisspeedbin))/sqrt(length(thisspeedbin));
    E06TISTDSPEED(I)=std(TI(thisspeedbin));
    E06TI90SPEED(I)=prctile(TI(thisspeedbin),90);
    clear  thisspeedbin
end
m_scatter_TI(E06_WS_160m_FLIDAR,E06_WS_TI_160m_FLIDAR,'ci',0:0.05:.25,'density','weights',ones(99999,1));
legend('Location','NorthEast')
hold on
errorbar(E06TISPEED,E06TISTDSPEED,'r','lineWidth',1,'DisplayName','TI')
plot([1:32],E06TI90SPEED,'r-','lineWidth',1,'DisplayName','TI 90')
[B,I]=sort(WS);
plot(WS(I),TI_offshore(I),'lineWidth',1,'DisplayName','TI offhsore');
plot(WS(I),TI_90_offshore(I),'lineWidth',1,'DisplayName','TI offshore 90');
plot(WS(I),TI_DNV_NTM_OA(I),'lineWidth',1,'DisplayName','DNV NTM OA');
plot(WS(I),TI_DNV_NTM_OB(I),'lineWidth',1,'DisplayName','DNV NTM OB');
plot(WS(I),TI_DNV_NTM_OC(I),'lineWidth',1,'DisplayName','DNV NTM OC');
set(gca,'Ylim',[0 .7]);
print -djpeg E06TIscatter.jpg;close
clear WS TI TI_offshore TI_90_offshore TI_IEC_NTM_A TI_IEC_NTM_B TI_IEC_NTM_C TI_DNV_NTM_OA TI_DNV_NTM_OB TI_DNV_NTM_OC
cd ..

%% E05 plot wave roses
cd ROSES
m_scatter(E05_waves_MWDHm0,E05_waves_Hm0,'di',[0.05 0.1 0.2 0.5:0.25:4.5])
%m_rose_plot(E05_waves_MWDHm0.data,E05_waves_Hm0.data,'Ay',[0:20:360],'di',[0.05 0.1 0.2 0.5:0.25:4.5])
%print -djpeg E05_wHm0_wave_rose.jpg;close
m_scatter(E05_waves_MWDHm0_swell,E05_waves_Hm0_swell,'di',[0.05 0.1 0.2 0.5:0.25:4.5])
%m_rose_plot(E05_waves_MWDHm0_swell.data,E05_waves_Hm0_swell.data,'Ay',[0:20:360],'di',[0.05 0.1 0.2 0.5:0.25:4.5])
%print -djpeg E05_wHm0_swell_wave_rose.jpg;close
m_scatter(E05_waves_MWDHm0_sea,E05_waves_Hm0_sea,'di',[0.05 0.1 0.2 0.5:0.25:4.5])
%m_rose_plot(E05_waves_MWDHm0_sea.data,E05_waves_Hm0_sea.data,'Ay',[0:20:360],'di',[0.05 0.1 0.2 0.5:0.25:4.5])
%print -djpeg E05_wHm0_sea_wave_rose.jpg;close

% E06 Wave Roses
% plot E06 wave roses
m_scatter(E06_waves_MWDHm0,E06_waves_Hm0,'di',[0.05 0.1 0.2 0.5:0.25:4.5])
m_scatter(E06_waves_MWDHm0_swell,E06_waves_Hm0_swell,'di',[0.05 0.1 0.2 0.5:0.25:4.5])
m_scatter(E06_waves_MWDHm0_sea,E06_waves_Hm0_sea,'di',[0.05 0.1 0.2 0.5:0.25:4.5])

% %m_rose_plot(E06_waves_MWDHm0.data,E06_waves_Hm0.data,'Ay',[0:20:360],'di',[0.05 0.1 0.2 0.5:0.25:4.5])
% %print -djpeg E06_wHm0_wave_rose.jpg;close
% m_rose_plot(E06_waves_MWDHm0_swell.data,E06_waves_Hm0_swell.data,'Ay',[0:20:360],'di',[0.05 0.1 0.2 0.5:0.25:4.5])
% print -djpeg E06_wHm0_swell_wave_rose.jpg;close
% m_rose_plot(E06_waves_MWDHm0_sea.data,E06_waves_Hm0_sea.data,'Ay',[0:20:360],'di',[0.05 0.1 0.2 0.5:0.25:4.5])
% print -djpeg E06_wHm0_sea_wave_rose.jpg;close
cd ..

%% E05 tidal analysis
cd TIDALANALYSIS
%E05_tidal_analysis=m_tide(E05_WL,'fmin',1/30);
E05_tidal_analysis=m_tide(E05_WL);
TidalLevels=m_tidal_levels(E05_tidal_analysis.WLTde)
save E05_tidal_analysis E05_tidal_analysis
% not sure what this does - called it and MAtlab was just hanging
%m_timeseries([E05_tidal_analysis.WL,E05_tidal_analysis.WLTde,E05_tidal_analysis.WLRes],'WriteTable',true);

% E06 tidal analysis
%E06_tidal_analysis=m_tide(E06_WL,'fmin',1/30);
E06_tidal_analysis=m_tide(E06_WL);
TidalLevels=m_tidal_levels(E06_tidal_analysis.WLTde)
save E06_tidal_analysis E06_tidal_analysis
E05tidedata=[E05_tidal_analysis.WL.data E05_tidal_analysis.WLTde.data E05_tidal_analysis.WLRes.data];
save E05tidedata.txt E05tidedata -ascii
E06tidedata=[E06_tidal_analysis.WL.data E06_tidal_analysis.WLTde.data E06_tidal_analysis.WLRes.data];
save E06tidedata.txt E06tidedata -ascii

cd ..

%% E05 plot current roses
cd ROSES
m_scatter(E05_current_direction,E05_current_speed,'di',[0 0.01 0.02 0.03 0.05:0.05:0.4]);
%m_rose_plot(E05_current_direction.data,E05_current_speed.data,'Ay',[0:20:360],'di',[0:0.05:0.4])
%print -djpeg E05_current_rose.jpg;close

%% E06 current  Roses
m_scatter(E06_current_direction,E06_current_speed,'di',[0 0.01 0.02 0.03 0.05:0.05:0.4]);
%m_rose_plot(E06_current_direction.data,E06_current_speed.data,'Ay',[0:20:360],'di',[0:0.05:0.4])
%print -djpeg E06_current_rose.jpg;close
cd ..

%% E05 current analysis
cd TIDALANALYSIS

%E05_current_analysis=m_tide(E05_current_speed,E05_current_direction,'fmin',1/30);
E05_current_analysis=m_tide(E05_current_speed,E05_current_direction);
fprintf(1,'%5.2f\n',prctile(E05_current_analysis.CS.data,[25 50 75 95])',max(E05_current_analysis.CS.data))
fprintf(1,'%5.2f\n',prctile(E05_current_analysis.CSRes.data,[25 50 75 95])',max(E05_current_analysis.CSRes.data))
fprintf(1,'%5.2f\n',prctile(E05_current_analysis.CSTde.data,[25 50 75 95])',max(E05_current_analysis.CSTde.data))
save E05_current_analysis E05_current_analysis
cd ..
cd ROSES
% raw data done above
m_scatter(E05_current_analysis.CDTde,E05_current_analysis.CSTde,'di',[0 0.01 0.02 0.03 0.05:0.05:0.4]);
m_scatter(E05_current_analysis.CDRes,E05_current_analysis.CSRes,'di',[0 0.01 0.02 0.03 0.05:0.05:0.4]);
% m_rose_plot(E05_current_analysis.CDTde.data,E05_current_analysis.CSTde.data,'Ay',[0:20:360],'di',[0.005 0.01 0.02 0.05:0.05:0.4])
% print -djpeg E05_current_rose_CSTde.jpg;close
% m_rose_plot(E05_current_analysis.CDRes.data,E05_current_analysis.CSRes.data,'Ay',[0:20:360],'di',[0.005 0.01 0.02 0.05:0.05:0.4])
% print -djpeg E05_current_rose_CSRes.jpg;close

E05currentdata=[E05_current_analysis.CS.data E05_current_analysis.CSTde.data E05_current_analysis.CSRes.data];
% export to txt
save E05currentdata.txt E05currentdata -ascii

%% E06 current analysis

%E06_current_analysis=m_tide(E06_current_speed,E06_current_direction,'fmin',1/30);
E06_current_analysis=m_tide(E06_current_speed,E06_current_direction);
fprintf(1,'%5.2f\n',prctile(E06_current_analysis.CS.data,[25 50 75 95])',max(E06_current_analysis.CS.data))
fprintf(1,'%5.2f\n',prctile(E06_current_analysis.CSRes.data,[25 50 75 95])',max(E06_current_analysis.CSRes.data))
fprintf(1,'%5.2f\n',prctile(E06_current_analysis.CSTde.data,[25 50 75 95])',max(E06_current_analysis.CSTde.data))
save E06_current_analysis E06_current_analysis 
cd ..
cd Roses
% raw data done above
m_scatter(E06_current_analysis.CDTde,E06_current_analysis.CSTde,'di',[0 0.01 0.02 0.03 0.05:0.05:0.4]);
m_scatter(E06_current_analysis.CDRes,E06_current_analysis.CSRes,'di',[0 0.01 0.02 0.03 0.05:0.05:0.4]);
%m_rose_plot(E06_current_analysis.CDTde.data,E06_current_analysis.CSTde.data,'Ay',[0:20:360],'di',[0.005 0.01 0.02 0.05:0.05:0.4])
%print -djpeg E06_current_rose_CSTde.jpg;close
%m_rose_plot(E06_current_analysis.CDRes.data,E06_current_analysis.CSRes.data,'Ay',[0:20:360],'di',[0.005 0.01 0.02 0.05:0.05:0.4])
%print -djpeg E06_current_rose_CSRes.jpg;close

E06currentdata=[E06_current_analysis.CS.data E06_current_analysis.CSTde.data E06_current_analysis.CSRes.data];
save E06currentdata.txt E06currentdata -ascii

cd ..

%% run this script to do wind profile plots and create tables of TI TI90
% bin by wind speed etc.
JWO_Wind_profile;

%% do air density
cd AIRDENSITY
E05_air_density=air_density_height(E05_temperature,E05_pressure,E05_RH,165);
for I=1:12
    thismonth=find(E05Airdate(:,2)==I);
    E05ADMONTH(I)=mean(E05_air_density(thismonth));
    E05ADSTDMONTH(I)=std(E05_air_density(thismonth))/sqrt(length(thismonth));
    clear thismonth
end
for I=1:29
    thisspeedbin=find(E05_CFSR_100m_speed.data>I-1 & E05_CFSR_100m_speed.data<I);
    E05ADSPEED(I)=mean(E05_air_density(thisspeedbin));
    E05ADSTDSPEED(I)=std(E05_air_density(thisspeedbin))/sqrt(length(thisspeedbin));
    clear  thisspeedbin
end


E06_air_density=air_density_height(E06_temperature,E06_pressure,E06_RH,165);
for I=1:12
    thismonth=find(E06Airdate(:,2)==I);
    E06ADMONTH(I)=mean(E06_air_density(thismonth));
    E06ADSTDMONTH(I)=std(E06_air_density(thismonth))/sqrt(length(thismonth));
    clear thismonth
end
for I=1:29
    thisspeedbin=find(E06_CFSR_100m_speed.data>I-1 & E06_CFSR_100m_speed.data<I);
    E06ADSPEED(I)=mean(E06_air_density(thisspeedbin));
    E06ADSTDSPEED(I)=std(E06_air_density(thisspeedbin))/sqrt(length(thisspeedbin));
    clear  thisspeedbin
end
% mean of all data
MEANAD=(mean(E05ADMONTH)+mean(E06ADMONTH))/2;
errorbar(E05ADMONTH,E05ADSTDMONTH,'b-')
hold on
errorbar(E06ADMONTH,E06ADSTDMONTH,'r-')
plot([1,12],[MEANAD,MEANAD],'k--');
set(gca,'xtick',1:12,'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
legend('E05','E06','Average')
grid on
xlabel('Month')
ylabel('Air Density (kg/m^3)')
print -djpeg Air_density_by_month.jpg;close

errorbar(E05ADSPEED,E05ADSTDSPEED,'b-')
hold on
errorbar(E06ADSPEED,E06ADSTDSPEED,'r-')
plot([1,35],[MEANAD,MEANAD],'k--');
legend('E05','E06','Average')
grid on
xlabel('Wind Speed (m/s)')
ylabel('Air Density (kg/m^3)')
print -djpeg Air_density_by_speed.jpg;close
cd ..

%% do scatter plots
cd SCATTER
% in the order in the report
% templates from report
% some bins need correcting based on Saraa's review comments - easieste to
% do it here rather than read in data again

% need to resample currents to one-hourly for scatter plots
CS_Total=jwo_resample(E05_current_speed,'ttt',[datenum([1979 01 01 01 00 00]) datenum([2020 12 31 23 00 00]) 30 60] ); 
% Total
CS_Total.bins=[0:0.1:1.0];
E05_waves_Hm0.bins=[0:1:15];
m_scatter(E05_waves_Hm0,CS_Total,'ci',0:0.05:.5,'quantiles',[0.05 0.5 0.95],'fit','Poly','density')
m_scatter(E05_waves_Hm0,E05_waves_Tp,'ci',0:5:15,'quantiles',[0.05 0.5 0.95],'fit','Poly','density')
E05_CFSR_165m_speed.label='U_{165}';
E05_CFSR_165m_speed.unit='m/s';
E05_CFSR_165m_speed.bins=[0:5:45];
m_scatter(E05_waves_Hm0,E05_CFSR_165m_speed,'ci',0:1:30,'quantiles',[0.05 0.5 0.95],'fit','Poly','density')
% Wind Sea
m_scatter(E05_waves_Hm0_sea,E05_CFSR_165m_speed,'ci',0:1:30,'quantiles',[0.05 0.5 0.95],'fit','Poly','density')
E05_waves_Hm0_sea.bins=[0:1:15];
m_scatter(E05_waves_Hm0_sea,CS_Total,'ci',0:0.05:.5,'quantiles',[0.05 0.5 0.95],'fit','Poly','density')
m_scatter(E05_waves_Hm0_sea,E05_waves_Tpsea,'ci',0:1:30,'quantiles',[0.05 0.5 0.95],'fit','Poly','density')
% Swell
E05_waves_Hm0_swell.bins=[0:1:15];
m_scatter(E05_waves_Hm0_swell,E05_CFSR_165m_speed,'ci',0:1:30,'quantiles',[0.05 0.5 0.95],'fit','Poly','density')
m_scatter(E05_waves_Hm0_swell,CS_Total,'ci',0:0.05:.5,'quantiles',[0.05 0.5 0.95],'fit','Poly','density')
m_scatter(E05_waves_Hm0_swell,E05_waves_Tpswell,'ci',0:1:30,'quantiles',[0.05 0.5 0.95],'fit','Poly','density')

% Air temp/windspeed
E05_CFSR10m_temperature.legend='CFSR 2 m';
E05_CFSR_10m_speed.bins=[0:2:30];
E05_CFSR_10m_temperature.bins=[-20:5:35];
m_scatter(E05_CFSR_10m_speed,E05_CFSR_10m_temperature,'quantiles',[0.05 0.5 0.95],'density')
E06_CFSR_10m_temperature.legend='CFSR 2 m';
E06_CFSR_10m_speed.bins=[0:2:30];
E06_CFSR_10m_temperature.bins=[-20:5:35];
m_scatter(E06_CFSR_10m_speed,E06_CFSR_10m_temperature,'quantiles',[0.05 0.5 0.95],'density')

clear CS_Total

% repeat for E06 as above 
CS_Total=jwo_resample(E06_current_speed,'ttt',[datenum([1979 01 01 01 00 00]) datenum([2020 12 31 23 00 00]) 30 60] ); 
% Total
CS_Total.bins=[0:0.1:1.0];
E06_waves_Hm0.bins=[0:1:15];
m_scatter(E06_waves_Hm0,CS_Total,'ci',0:0.05:.5,'quantiles',[0.05 0.5 0.95],'fit','Poly','density')
m_scatter(E06_waves_Hm0,E06_waves_Tp,'ci',0:5:15,'quantiles',[0.05 0.5 0.95],'fit','Poly','density')
E06_CFSR_165m_speed.label='U_{165}';
E06_CFSR_165m_speed.unit='m/s';
E06_CFSR_165m_speed.bins=[0:5:45];
m_scatter(E06_waves_Hm0,E06_CFSR_165m_speed,'ci',0:1:30,'quantiles',[0.05 0.5 0.95],'fit','Poly','density')
% Wind Sea
m_scatter(E06_waves_Hm0_sea,E06_CFSR_165m_speed,'ci',0:1:30,'quantiles',[0.05 0.5 0.95],'fit','Poly','density')
E06_waves_Hm0_sea.bins=[0:1:15];
m_scatter(E06_waves_Hm0_sea,CS_Total,'ci',0:0.05:.5,'quantiles',[0.05 0.5 0.95],'fit','Poly','density')
m_scatter(E06_waves_Hm0_sea,E06_waves_Tpsea,'ci',0:1:30,'quantiles',[0.05 0.5 0.95],'fit','Poly','density')
% Swell
E06_waves_Hm0_swell.bins=[0:1:15];
m_scatter(E06_waves_Hm0_swell,E06_CFSR_165m_speed,'ci',0:1:30,'quantiles',[0.05 0.5 0.95],'fit','Poly','density')
m_scatter(E06_waves_Hm0_swell,CS_Total,'ci',0:0.05:.5,'quantiles',[0.05 0.5 0.95],'fit','Poly','density')
m_scatter(E06_waves_Hm0_swell,E06_waves_Tpswell,'ci',0:1:30,'quantiles',[0.05 0.5 0.95],'fit','Poly','density')


clear CS_Total

cd ..


%% misalignment - not sure how to do density and also seperate plot for < 20 and 20+
cd MISALIGNMENT
% waves and 10m and 165m winds
m_misalignment(E05_waves_MWDHm0,E05_CFSR_10m_direction,E05_CFSR_10m_speed,'directional',E05_CFSR_10m_direction);
m_misalignment(E06_waves_MWDHm0,E06_CFSR_10m_direction,E06_CFSR_10m_speed,'directional',E06_CFSR_10m_direction);
m_misalignment(E05_waves_MWDHm0,E05_CFSR_165m_direction,E05_CFSR_165m_speed,'directional',E05_CFSR_165m_direction);
m_misalignment(E06_waves_MWDHm0,E06_CFSR_165m_direction,E06_CFSR_165m_speed,'directional',E06_CFSR_165m_direction);
% currents and 10m and 165m wind
CD_Total=jwo_resample(E05_current_direction,'ttt',[datenum([1979 01 01 01 00 00]) datenum([2020 12 31 23 00 00]) 30 60]); 
m_misalignment(CD_Total,E05_CFSR_10m_direction,E05_CFSR_10m_speed,'directional',CD_Total);
m_misalignment(CD_Total,E05_CFSR_165m_direction,E05_CFSR_165m_speed,'directional',CD_Total);
clear CD_Total
CD_Total=jwo_resample(E06_current_direction,'ttt',[datenum([1979 01 01 01 00 00]) datenum([2020 12 31 23 00 00]) 30 60]); 
m_misalignment(CD_Total,E06_CFSR_10m_direction,E06_CFSR_10m_speed,'directional',CD_Total);
m_misalignment(CD_Total,E06_CFSR_165m_direction,E06_CFSR_165m_speed,'directional',CD_Total);
clear CD_Total
cd ..

%% swell assessment plots
cd SWELL
% in the order in the report
m_scatter(E05_waves_Hm0_sea,E05_waves_Hm0,'ci',0:5:15,'quantiles',[0.05 0.5 0.95],'fit','Poly','density')
m_scatter(E05_waves_Hm0_swell,E05_waves_Hm0_sea,'ci',0:5:15,'quantiles',[0.05 0.5 0.95],'fit','Poly','density')
m_scatter(E05_waves_Hm0_swell,E05_waves_Hm0,'ci',0:5:15,'quantiles',[0.05 0.5 0.95],'fit','Poly','density')

m_scatter(E06_waves_Hm0_sea,E06_waves_Hm0,'ci',0:5:15,'quantiles',[0.05 0.5 0.95],'fit','Poly','density')
m_scatter(E06_waves_Hm0_swell,E06_waves_Hm0_sea,'ci',0:5:15,'quantiles',[0.05 0.5 0.95],'fit','Poly','density')
m_scatter(E06_waves_Hm0_swell,E06_waves_Hm0,'ci',0:5:15,'quantiles',[0.05 0.5 0.95],'fit','Poly','density')
cd ..

%% Persistence
% save the variables neeeded - so that a seperate script can be run 
% cycling through parameters, percentile, thresholds, durations
cd PERSISTENCE
save PersistenceData E05_CFSR_10m_speed E06_CFSR_10m_speed E05_waves_Hm0 E06_waves_Hm0 E05_waves_Tp E06_waves_Tp 
% example for report
HDUR_WIND  = [1 6 12 18 24 30 36 42 48 54 60 66 72]; % h
THR_WIND  = [1 3 5 7 9 11 13 15 17 19]; % m/s
m_persistence(E05_CFSR_10m_speed,'DUR',HDUR_WIND,'THR',THR_WIND,'P',50,'plot_flag',[1]);
cd ..

