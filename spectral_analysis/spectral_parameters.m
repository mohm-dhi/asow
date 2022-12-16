
clear all; close all; clc

load('\\USDEN1-STOR.DHI.DK\Projects\41806529\02_RAW_MOOD_data\2D Spectra data\Spec0.1_1979-2021_-73.9W39.3N.mat','X')

X.dfs_unit = X.unit;
X.dfs_info.UnitAbbr = 'm^2*s/deg';

%%
%X.time = X.time(1:1000);
%X.data = X.data(1:1000,:,:);

IntParam = m_2DspecInt(X);
% IntParam = [Hm0,Tp,T01,T02,Tm10,PWD,MWD,DSD,DSD_p,Spr,Spr_p,PrWD,Tp2,Spr_p2,PWD2];

params_fname = '\\USDEN1-STOR.DHI.DK\Projects\41806529\07_Timeseries_CSV_Deliverables\ListParameters.xlsx';
translateItems = readtable(params_fname,'Sheet','TranslateSpectralData');

S = struct();
S.Time = datetime(IntParam.Time,'ConvertFrom','datenum');
for k = 1:height(translateItems)
    
    mask = strcmp(IntParam.ItemNames,translateItems.Spec{k});
    
    S.(translateItems.Final{k}) = IntParam.Values(:,mask);
    
end

S.xyz = [-73.9,39.3,33.4];

save('\\USDEN1-STOR.DHI.DK\Projects\41806529\03_PROCESSED_MOOD_data\02-OutputLocations\SpectralAnalysis\Spectra_Integral_Parameters_Loc008.mat','IntParam','S');

%%

figure('Color','w','POsition',[0 0 800 500]); 
subplot(121);hold on
cmap = hot(30);
cmap = flipud(cmap(5:end-1,:));
colormap(cmap)
%plot(S.Tp_Total, S.T10_Total,'k.','markersize',2);
hold on
histogram2(S.Tp_Total, S.T10_Total,0:25,0:25,'DisplayStyle','tile','EdgeColor','None');
plot([2, 20],[2, 20]*0.81,'DisplayName','Tm-1,0 = 0.81Tp')
plot([2, 20],[2, 20]*0.97,'DisplayName','Tm-1,0 = 0.97Tp')

daspect([1 1 1])
%legend()
xlabel('Tp [s]')
ylabel('Tm-1,0 [s]')
grid on; box on;
subplot(122);hold on
cmap = hot(30);
cmap = flipud(cmap(5:end-1,:));
colormap(cmap)
plot(S.Tp_Total, S.T10_Total,'.','markersize',2,'HandleVisibility','off');
hold on
%histogram2(S.Tp_Total, S.T10_Total,0:25,0:25,'DisplayStyle','tile','EdgeColor','None');
plot([2, 20],[2, 20]*0.81,'DisplayName','Tm-1,0 = 0.81Tp')
plot([2, 20],[2, 20]*0.97,'DisplayName','Tm-1,0 = 0.97Tp')

daspect([1 1 1])
axis([0 25 0 25])
legend('Location','NW')
xlabel('Tp [s]')
ylabel('Tm-1,0 [s]')
grid on; box on;

print('Comparison_Periods.png','-dpng','-r300')