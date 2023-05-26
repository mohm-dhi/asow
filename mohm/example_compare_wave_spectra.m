% Example on plotting wave spectra model and measurements for comparison
% Adopted from N-7.2 OWF project (11827978)
% Date: 2023-01-16
% Author: ARR
%

clc
clear all
close all

% options
opt.print_fig = false;

% load mtools in path
addpath(genpath('\\dkcph1-stor\projects\11827978\_COMMON_DATA_\MATLAB\potlab_v2_copy\res\DHI-MATLAB-Toolbox'))
addpath(genpath('\\dkcph1-stor\projects\11827978\_COMMON_DATA_\MATLAB\potlab_v2_copy\res\wafo'))
addpath(genpath('\\dkcph1-stor\projects\11827978\_COMMON_DATA_\MATLAB\potlab_v2_copy\res\UTide'))
addpath(genpath('\\dkcph1-stor\projects\11827978\_COMMON_DATA_\MATLAB\potlab_v2_copy\src\m_tools'))

% find dfs2 wave spectra (modelled)
rf='\\DKCPH1-STOR2.DHI.DK\Projects\11827978_WS\_COMMON_\recover_dwf2020\swdwf2020_spectra_5km_n72\';
dfs2_list=dir([rf '*.dfs2']);

%% spectra comparison plot
ttt=[datenum([2017 9 16]) datenum([2020 6 1]) 60];
%load model spectra
model_spec_fn='\\dkcph1-stor\projects\11827978\CHC\Model Validation\Offshore_1976_f-spec.dfs1';
% write into m_structure
SW_Spec      = m_structure('Offshore 1976',[6.33,54.28,-30],ttt ,'SW_{DWF2020}',model_spec_fn      ,'ED1f'     , 1,0:10:200);

%load measured spectra
% NOTE: dfs1 dont accept non-equidistant because of gap in measurement so split into 3 files
meas_spec1_fn='\\dkcph1-stor\projects\11827978\CHC\Measurements\N7S\spectra_1.dfs1';
meas_spec2_fn='\\dkcph1-stor\projects\11827978\CHC\Measurements\N7S\spectra_2.dfs1';
meas_spec3_fn='\\dkcph1-stor\projects\11827978\CHC\Measurements\N7S\spectra_3.dfs1';

meas_Spec1      = m_structure('Offshore 1976',[6.33,54.28,-30],ttt ,'SW_{DWF2020}',meas_spec1_fn      ,'ED1f'     , 1,0:10:200);
meas_Spec2      = m_structure('Offshore 1976',[6.33,54.28,-30],ttt ,'SW_{DWF2020}',meas_spec2_fn      ,'ED1f'     , 1,0:10:200);
meas_Spec3      = m_structure('Offshore 1976',[6.33,54.28,-30],ttt ,'SW_{DWF2020}',meas_spec3_fn      ,'ED1f'     , 1,0:10:200);

%merge all the measured spectra
meas_Spec_all=meas_Spec1;
meas_Spec_all.data=[meas_Spec1.data;meas_Spec2.data;meas_Spec3.data];
meas_Spec_all.time=[meas_Spec1.time;meas_Spec2.time;meas_Spec3.time];

MeSavg1=nanmean(meas_Spec1.data);
MeSavg2=nanmean(meas_Spec2.data);
MeSavg3=nanmean(meas_Spec3.data);

%plot density comparison
FigH = figure('Position', get(0, 'Screensize'));
subplot(3,2,1);
[X1,Y1]=meshgrid(meas_Spec1.time,meas_Spec1.xaxis);
pcolor(X1',Y1',meas_Spec1.data)
shading interp
caxis([0 20])
colorbar('FontSize',12)
xlabel('Date (yyyy/mm/dd)','FontSize',12)
ylabel('Frequency (Hz)','FontSize',12)
axis([meas_Spec1.time(1) meas_Spec1.time(end) 0 0.4])
datetick('x','yyyy/mm/dd','keeplimits')
title('spectral energy (measurement) [m^{2}/Hz]','FontSize',12)

subplot(3,2,2);
[X1,Y1]=meshgrid(SW_Spec.time,SW_Spec.xaxis);
pcolor(X1',Y1',SW_Spec.data)
shading interp
caxis([0 20])
colorbar('FontSize',12)
xlabel('Date (yyyy/mm/dd)','FontSize',12)
ylabel('Frequency (Hz)','FontSize',12)
axis([meas_Spec1.time(1) meas_Spec1.time(end) 0 0.4])
datetick('x','yyyy/mm/dd','keeplimits')
title('spectral energy (model) [m^{2}/Hz]','FontSize',12)

subplot(3,2,3);
[X1,Y1]=meshgrid(meas_Spec2.time,meas_Spec2.xaxis);
pcolor(X1',Y1',meas_Spec2.data)
shading interp
caxis([0 20])
colorbar('FontSize',12)
xlabel('Date (yyyy/mm/dd)','FontSize',12)
ylabel('Frequency (Hz)','FontSize',12)
axis([meas_Spec2.time(1) meas_Spec2.time(end) 0 0.4])
datetick('x','yyyy/mm/dd','keeplimits')
title('spectral energy (measurement) [m^{2}/Hz]','FontSize',12)

subplot(3,2,4);
[X1,Y1]=meshgrid(SW_Spec.time,SW_Spec.xaxis);
pcolor(X1',Y1',SW_Spec.data)
shading interp
caxis([0 20])
colorbar('FontSize',12)
xlabel('Date (yyyy/mm/dd)','FontSize',12)
ylabel('Frequency (Hz)','FontSize',12)
axis([meas_Spec2.time(1) meas_Spec2.time(end) 0 0.4])
datetick('x','yyyy/mm/dd','keeplimits')
title('spectral energy (model) [m^{2}/Hz]','FontSize',12)

subplot(3,2,5);
[X1,Y1]=meshgrid(meas_Spec3.time,meas_Spec3.xaxis);
pcolor(X1',Y1',meas_Spec3.data)
shading interp
caxis([0 20])
colorbar('FontSize',12)
xlabel('Date (yyyy/mm/dd)','FontSize',12)
ylabel('Frequency (Hz)','FontSize',12)
axis([meas_Spec3.time(1) meas_Spec3.time(end) 0 0.4])
datetick('x','yyyy/mm/dd','keeplimits')
title('spectral energy (measurement) [m^{2}/Hz]','FontSize',12)

subplot(3,2,6);
[X1,Y1]=meshgrid(SW_Spec.time,SW_Spec.xaxis);
pcolor(X1',Y1',SW_Spec.data)
shading interp
caxis([0 20])
colorbar('FontSize',12)
xlabel('Date (yyyy/mm/dd)','FontSize',12)
ylabel('Frequency (Hz)','FontSize',12)
axis([meas_Spec3.time(1) meas_Spec3.time(end) 0 0.4])
datetick('x','yyyy/mm/dd','keeplimits')
title('spectral energy (model) [m^{2}/Hz]','FontSize',12)

if opt.print_fig == true
    print(FigH,'spectra_comparison.png','-dpng','-r300');
end

[c1sta i1sta] = min(abs(SW_Spec.time-meas_Spec1.time(1)));
[c1end i1end] = min(abs(SW_Spec.time-meas_Spec1.time(end)));
[c2sta i2sta] = min(abs(SW_Spec.time-meas_Spec2.time(1)));
[c2end i2end] = min(abs(SW_Spec.time-meas_Spec2.time(end)));
[c3sta i3sta] = min(abs(SW_Spec.time-meas_Spec3.time(1)));
[c3end i3end] = min(abs(SW_Spec.time-meas_Spec3.time(end)));

MoSavg1=mean(SW_Spec.data(i1sta:i1end,:));
MoSavg2=mean(SW_Spec.data(i2sta:i2end,:));
MoSavg3=mean(SW_Spec.data(i3sta:i3end,:));

FigHX = figure('Position', get(0, 'Screensize'));
subplot(3,1,1);
plot(meas_Spec1.xaxis,MeSavg1,'LineWidth',3); hold on;
plot(SW_Spec.xaxis,MoSavg1,'LineWidth',3); 
xlabel('frequency [Hz]','FontSize',12)
ylabel('spectral energy [m^{2}/Hz]','FontSize',12)
legend('measurement','model','FontSize',12)
title('averaged spectral energy [m^{2}/Hz] (2017/9/16 to 2018/2/3)','FontSize',12)    
subplot(3,1,2);
plot(meas_Spec2.xaxis,MeSavg2,'LineWidth',3); hold on;
plot(SW_Spec.xaxis,MoSavg2,'LineWidth',3); 
xlabel('frequency [Hz]','FontSize',12)
ylabel('spectral energy [m^{2}/Hz]','FontSize',12)
legend('measurement','model','FontSize',12)
title('averaged spectral energy [m^{2}/Hz] (2018/7/31 to 2019/5/8)','FontSize',12)   
subplot(3,1,3);
plot(meas_Spec3.xaxis,MeSavg3,'LineWidth',3); hold on;
plot(SW_Spec.xaxis,MoSavg3,'LineWidth',3); 
xlabel('frequency [Hz]','FontSize',12)
ylabel('spectral energy [m^{2}/Hz]','FontSize',12)
legend('measurement','model','FontSize',12)
title('averaged spectral energy [m^{2}/Hz] (2020/1/21  to 2020/5/31)','FontSize',12)  

if opt.print_fig == true
    %save figure
    print(FigHX,'SpectraCompareAvg3.png','-dpng','-r300');
end