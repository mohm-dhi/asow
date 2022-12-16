% check on CNT and H0 TI_model script - used values from m_scatter_turb.m
% load in the apprpropriate windspeed and standard deviation data first
% from csv file of 10 minute FLIDAR data at 160m with all NaN's removed
clc
clear all
% read in speed, std, direction and reported TI from FLIDAR
fname='E05_Hudson_North_10_min_avg_20190812_20210919_160m.csv';
E05_WS_160m_FLIDAR = m_structure('E05' , [-72.72 39.97 57] , [datenum(2019,08,12,0,0,0) datenum(2021,9,19,23,50,0) 10] , 'FLIDAR Speed' , fname ,'U160' ,1,[0:2:34]);
E05_WS_Std_160m_FLIDAR = m_structure('E05' , [-72.72 39.97 57] , [datenum(2019,08,12,0,0,0) datenum(2021,9,19,23,50,0) 10] , 'FLIDAR Std' , fname ,'WSSTD160' ,2,[0:.1:2]);
E05_WD_160m_FLIDAR = m_structure('E05' , [-72.72 39.97 57] , [datenum(2019,08,12,0,0,0) datenum(2021,9,19,23,50,0) 10] , 'FLIDAR Direction' , fname ,'D160' ,3,[0:30:360]);
E05_WS_TI_160m_FLIDAR = m_structure('E05' , [-72.72 39.97 57] , [datenum(2019,08,12,0,0,0) datenum(2021,9,19,23,50,0) 10] , 'FLIDAR TI' , fname ,'TI160',4,[0:.1:2]);
% replace TI in raw data with std/V ??
E05_WS_TI_160m_FLIDAR.data=E05_WS_Std_160m_FLIDAR.data./E05_WS_160m_FLIDAR.data;

fname='E06_Hudson_South_10_min_avg_20190904_20220123_160m.csv';
E06_WS_160m_FLIDAR = m_structure('E06' , [-72.72 39.97 34] , [datenum(2019,9,4,0,0,0) datenum(2022,1,23,23,50,0) 10] , 'FLIDAR WS' , fname ,'U160' ,1,[0:2:34]);
E06_WS_Std_160m_FLIDAR = m_structure('E06' , [-72.72 39.97 34] , [datenum(2019,9,4,0,0,0) datenum(2022,1,23,23,50,0) 10] , 'FLIDAR Std' , fname ,'WSSTD160' ,2,[0:.1:2]);
E06_WD_160m_FLIDAR = m_structure('E06' , [-72.72 39.97 34] , [datenum(2019,9,4,0,0,0) datenum(2022,1,23,23,50,0) 10] , 'FLIDAR Direction' , fname ,'D160' ,3,[0:30:360]);
E06_WS_TI_160m_FLIDAR = m_structure('E06' , [-72.72 39.97 34] , [datenum(2019,9,4,0,0,0) datenum(2022,1,23,23,50,0) 10] , 'FLIDAR TI' , fname ,'TI160' ,4,[0:.1:2]);
% replace TI in raw data with std/V ??
E06_WS_TI_160m_FLIDAR.data=E06_WS_Std_160m_FLIDAR.data./E06_WS_160m_FLIDAR.data;
% % messy but needs m_structure call for m_scatter to work!
% 
save JWOTI

load JWOTI

WS=E05_WS_160m_FLIDAR.data';
TI=E05_WS_TI_160m_FLIDAR.data';

m_rose_plot(E05_WD_160m_FLIDAR.data,E05_WS_TI_160m_FLIDAR.data,'Ay',[0:30:360],'ci',0:5:15,'di',[0 0.01 0.02 0.05 0.1 0.2 0.3],'labtitle',['E05 TI at 160m'] )
print -djpeg E05TIrose.jpg;close

% call TI_models
[TI_offshore,TI_90_offshore,TI_IEC_NTM_A,TI_IEC_NTM_B,TI_IEC_NTM_C,TI_DNV_NTM_OA,TI_DNV_NTM_OB,TI_DNV_NTM_OC] = TI_models(165,WS,2.4,24.2);

% split into month and time of day
for I=1:34
    thisspeedbin=find(WS>I-1 & WS<I);
    E05TISPEED(I)=mean(TI(thisspeedbin));
    E05TISTDSPEED(I)=std(TI(thisspeedbin));
    E05TISTDSPEED(I)=std(TI(thisspeedbin))/sqrt(length(thisspeedbin));
    E05TI90SPEED(I)=prctile(TI(thisspeedbin),90);
    clear  thisspeedbin
end

% do coloued scatter plot
m_scatter_TI(E05_WS_160m_FLIDAR,E05_WS_TI_160m_FLIDAR,'ci',0:0.05:.25,'density','weights',ones(99999,1));
legend('Location','NorthEast')
hold on
% workaround to plot TI estimates
%errorbar(E05TISPEED,E05TISTDSPEED,'r','lineWidth',1,'DisplayName','TI')
plot([1:34],E05TISPEED,'r-','lineWidth',1,'DisplayName','TI')
plot([1:34],E05TI90SPEED,'r--','lineWidth',1,'DisplayName','TI 90')
% sort the estiamtes into order so plot looks like a single line
[B,I]=sort(WS);
plot(WS(I),TI_offshore(I),'lineWidth',1,'DisplayName','TI offshore');
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
for I=1:34
    thisspeedbin=find(WS>I-1 & WS<I);
    E06TISPEED(I)=mean(TI(thisspeedbin));
    E06TISTDSPEED(I)=std(TI(thisspeedbin));
    E06TISTDSPEED(I)=std(TI(thisspeedbin))/sqrt(length(thisspeedbin));
    E06TI90SPEED(I)=prctile(TI(thisspeedbin),90);
    clear  thisspeedbin
end
m_scatter_TI(E06_WS_160m_FLIDAR,E06_WS_TI_160m_FLIDAR,'ci',0:0.05:.25,'density','weights',ones(99999,1));
legend('Location','NorthEast')
hold on
%errorbar(E06TISPEED,E06TISTDSPEED,'r','lineWidth',1,'DisplayName','TI')
plot([1:34],E06TISPEED,'r-','lineWidth',1,'DisplayName','TI')
plot([1:34],E06TI90SPEED,'r--','lineWidth',1,'DisplayName','TI 90')
[B,I]=sort(WS);
plot(WS(I),TI_offshore(I),'lineWidth',1,'DisplayName','TI offshore');
plot(WS(I),TI_90_offshore(I),'lineWidth',1,'DisplayName','TI offshore 90');
plot(WS(I),TI_DNV_NTM_OA(I),'lineWidth',1,'DisplayName','DNV NTM OA');
plot(WS(I),TI_DNV_NTM_OB(I),'lineWidth',1,'DisplayName','DNV NTM OB');
plot(WS(I),TI_DNV_NTM_OC(I),'lineWidth',1,'DisplayName','DNV NTM OC');
set(gca,'Ylim',[0 .7]);
print -djpeg E06TIscatter.jpg;close
clear WS TI TI_offshore TI_90_offshore TI_IEC_NTM_A TI_IEC_NTM_B TI_IEC_NTM_C TI_DNV_NTM_OA TI_DNV_NTM_OB TI_DNV_NTM_OC