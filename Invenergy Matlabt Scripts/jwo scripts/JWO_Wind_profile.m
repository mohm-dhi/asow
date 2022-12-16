% read in all speed and std data from 60 to 200 m
% plot TI/TI90 vs WS
fname= '..\NYSERDAFLiDAR\Data\E05_Hudson_North_10_min_avg_20190812_20210919.csv';
RawFlidar  = table2array(readtable(fname));
E05WSRAW=RawFlidar(:,[7:10:97]);
E05STDRAW=RawFlidar(:,[8:10:98]);
E05WDRAW=RawFlidar(:,[11:10:101]);
E05TI=E05STDRAW./E05WSRAW;

%TI rose plots
m_rose_plot(E05WDRAW(:,8),E05TI(:,8),'Ay',[0:30:360],'di',[0.0:0.03:0.3])
print -djpeg E05_TI_rose.jpg;close
% bin by wind speed
for I=1:28
    SpeedBin=find(E05WSRAW(:,8)>=I & E05WSRAW(:,8)<I+1);
    E05TITABLE(I,1)=mean(E05TI(SpeedBin,8));
    E05TITABLE(I,2)=prctile(E05TI(SpeedBin,8),90);
    E05TITABLE(I,3)=std(E05TI(SpeedBin,8));
    clear SpeedBin
end

% extract mean TI/TI90 from 60-200 m
for J=3:10
    for I=1:35
        SpeedBin=find(E05WSRAW(:,J)>=I & E05WSRAW(:,J)<I+1);
        E05MeanTI(I,J-2)=mean(E05TI(SpeedBin,J));
        E05StdTI(I,J-2)=std(E05TI(SpeedBin,J))/sqrt(length(SpeedBin));
        E05TI90(I,J-2)=prctile(E05TI(SpeedBin,J),90);
        clear SpeedBin
    end
end
% do the plot
for I=1:8
    errorbar([1:35],E05MeanTI(:,I),E05StdTI(:,I))
    hold on
end
grid on
xlabel('WS (m/s)')
ylabel('TI[-]')
title('E05 TI by height (binned by windspeed)')
set(gca,'Ylim',[0 0.25])
legend('60m','80m','100,','120m','140m','160m','180m','200m')
print -djpeg E05_TI_byheight.jpg;close

figure
plot(E05TI90)
grid on
xlabel('WS (m/s)')
ylabel('TI_9_0[-]')
set(gca,'Ylim',[0 0.25])
legend('60m','80m','100,','120m','140m','160m','180m','200m')
title('E06 TI_9_0 by height (binned by windspeed)')
print -djpeg E05_TI90_byheight.jpg;close

clear RawFlidar


fname= '..\NYSERDAFLiDAR\Data\E06_Hudson_South_10_min_avg_20190904_20220123.csv';
RawFlidar  = table2array(readtable(fname));
E06WSRAW=RawFlidar(:,[7:10:97]);
E06STDRAW=RawFlidar(:,[8:10:98]);
E06WDRAW=RawFlidar(:,[11:10:101]);
E06TI=E06STDRAW./E06WSRAW;

%TI rose plots
m_rose_plot(E06WDRAW(:,8),E06TI(:,8),'Ay',[0:30:360],'di',[0.0:0.03:0.3])
print -djpeg E06_TI_rose.jpg;close
% bin by wind speed
for I=1:28
    SpeedBin=find(E06WSRAW(:,8)>=I & E06WSRAW(:,8)<I+1);
    E06TITABLE(I,1)=mean(E06TI(SpeedBin,8));
    E06TITABLE(I,2)=prctile(E06TI(SpeedBin,8),90);
    E06TITABLE(I,3)=std(E06TI(SpeedBin,8));
    clear SpeedBin
end

% extract mean TI/TI90 from 60-200 m
for J=3:10
    for I=1:35
        SpeedBin=find(E06WSRAW(:,J)>=I & E06WSRAW(:,J)<I+1);
        E06MeanTI(I,J-2)=mean(E06TI(SpeedBin,J));
        E06StdTI(I,J-2)=std(E06TI(SpeedBin,J))/sqrt(length(SpeedBin));
        E06TI90(I,J-2)=prctile(E06TI(SpeedBin,J),90);
        clear SpeedBin
    end
end
%do the plot
for I=1:8
    errorbar([1:35],E06MeanTI(:,I),E06StdTI(:,I))
    hold on
end
grid on
xlabel('WS (m/s)')
ylabel('TI[-]')
set(gca,'Ylim',[0 0.25])
legend('60m','80m','100,','120m','140m','160m','180m','200m')
title('E06 TI by height (binned by windspeed)')
print -djpeg E06_TI_byheight.jpg;close

figure
plot(E06TI90)
grid on
xlabel('WS (m/s)')
ylabel('TI_9_0[-]')
set(gca,'Ylim',[0 0.25])
legend('60m','80m','100,','120m','140m','160m','180m','200m')
title('E06 TI_9_0 by height (binned by windspeed)')
print -djpeg E06_TI90_byheight.jpg;close

subplot(211);scatter(E05WSRAW(:,8),E05TI(:,8),'k.');axis square
subplot(212);scatter(E06WSRAW(:,8),E06TI(:,8),'k.');axis square
print -djpeg TI_scatter.jpg;close

clear I J fname RawFlidar


% plot mean profile
subplot(211)
loglog(nanmean(E05WSRAW),[20:20:200],'k*-')
set(gca,'Xlim',[8 12])
set(gca,'Ylim',[20 180])
grid on
xlabel('Mean windspeed (m/s) at site E05')
ylabel('Height (m)')
grid on
subplot(212)
loglog(nanmean(E06WSRAW),[20:20:200],'k*-')
set(gca,'Xlim',[8 12])
set(gca,'Ylim',[20 180])
grid on
xlabel('Mean windspeed (m/s) at site E06')
ylabel('Height (m)')
grid on

print -djpeg meanwindprofiles.jpg;close
% rose plots
for I=1:10
    m_rose_plot(E05WDRAW(:,I),E05WSRAW(:,I),'Ay',[0:30:360],'di',[0:2:24],'labtitle',['E05 FLIDAR at ' num2str(20+(I-1)*20) ' m'] )
    eval(['print -djpeg E05' num2str(20+(I-1)*20) 'm.jpg;close'])
end
for I=1:10
    m_rose_plot(E06WDRAW(:,I),E06WSRAW(:,I),'Ay',[0:30:360],'di',[0:2:24],'labtitle',['E06 FLIDAR at ' num2str(20+(I-1)*20) ' m'] )
    eval(['print -djpeg E06' num2str(20+(I-1)*20) 'm.jpg;close'])
end