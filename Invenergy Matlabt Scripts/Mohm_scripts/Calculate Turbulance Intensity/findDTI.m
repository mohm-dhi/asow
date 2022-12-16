close all;
clear all;
% %find delta sigma for TI
% figure
% V=output{1}{3}.X
% sigma=(output{1}{3}.Q2-output{1}{3}.Q1).*V
% plot(V,sigma,'.r')
% hold on
% fit=polyfit(V,sigma,1)
% plot(V,fit(1)*V+fit(2),'-b')
% hold off
%
% Vhub=[0:1:50];

sites = {'E05', 'E06'};
file_paths = {
    'C:\Users\mohm\OneDrive - DHI\Documents\Project Files\41806038 - Invenergy_NY Bight_Metocean\NYBightData\CFSR Hub Height\CFSR Hub 165 E05\CFSR_at_hubHeight_165m.dfs0'
    'C:\Users\mohm\OneDrive - DHI\Documents\Project Files\41806038 - Invenergy_NY Bight_Metocean\NYBightData\CFSR Hub Height\CFSR Hub 165 E06\CFSR_at_hubHeight_165m.dfs0'
    };

TI_files = {
    'C:\Users\mohm\OneDrive - DHI\Documents\Project Files\41806038 - Invenergy_NY Bight_Metocean\NYBightData\NYSERDAFLiDAR\Data\E05_Hudson_North_10_min_avg_20190812_20210919.csv'
    'C:\Users\mohm\OneDrive - DHI\Documents\Project Files\41806038 - Invenergy_NY Bight_Metocean\NYBightData\NYSERDAFLiDAR\Data\E06_Hudson_South_10_min_avg_20190904_20220123.csv';
    
    };
%% read data
for f = 1: length(file_paths)
    site = sites{f};
    file_path = file_paths{f};
    TI_file = TI_files{f};
    [~,data] = dfs02mat_dotnet(file_path);
    TI_obs = xlsread(TI_file);
    %% clean data
    
    TI = TI_obs(:,80); % data at near hub height 160m
    WS = TI_obs(:,71); % data at near hub height 160m
    
    I = ~isnan(TI) & ~isnan(WS);
    TI = TI(I);
    WS = WS(I);
    
    
    fig=figure;
    Vhub = sort(data.Values(:,1));
    Vave=mean(Vhub);
    Vave = 10.35 ;
    Iref = 0.1;
    Iref_dnv = 0.72;
    sigmaETM1=1.4*Iref*((3.*Vhub+38)/4.-(Vhub-Vave)/18.);
    
    sigmaETM2_ti=TI+4.5.*std(TI);
    
    ff=12
    plot(Vhub,sigmaETM1./Vhub,'r.')
    hold on
    ws = 0:29;
    for i=1:length(ws)-1
        
        sigma_ti(i) = mean(TI(WS>ws(i) & WS<= ws(i+1))) + 4.5* std(TI(WS>ws(i) & WS<= ws(i+1)));
        
    end
    
    
    plot(1:29,sigma_ti,'.b','MarkerSize',14)
    plot(Vhub,Iref_dnv*sigmaETM1./Vhub,'--b')
    ylim([0,0.5])
    title(['Offshore extreme turbulence model at ' site],'FontSize',ff)
    ylabel('TI_{OETM} [-]','FontSize',ff)
    xlabel('wind speed [m/s]','FontSize',ff)
    legend(['DNV standard I_{ref}=' num2str(Iref)],'ASIT LiDAR data',['DNV standard I_{ref}=' num2str(Iref*Iref_dnv) ' (fit)'],'FontSize',ff)
    
    saveas(fig,['oetm' site '.png'])
    
end
