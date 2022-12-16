clc; 
close all;
clear all;
warning off;

addpath(genpath('C:\Users\JABE\GitHub\potlab\third'))
addpath(genpath('C:\Users\JABE\GitHub\potlab\toolboxes\m_tools'))
addpath(genpath('C:\Users\JABE\GitHub\potlab\users\jabe'))

load('MetMast.mat')
[weight_months,weight_days] = monthly_daily_means(T.Time,6);
T=addvars(T,weight_months,weight_days);



size_crit=18;

fs=13;

heights=[25,50,70,80,90,100];
boom_dir='NW';
boom_dir1='NW';
boom_dir2='NW';

print_label='ABS_NW';

abs_on=1;

c=1;
zref=5;
var_dir_ref=strcat(boom_dir,'_WD_100m_Avg');
var_speed_ref=strcat(boom_dir,'_WS_',num2str(heights(zref)),'m_Avg');;
z1=1; 
z2=6;

hFit_min=1;
hFit_max=z2-z1+1;
 

 
nSectors=12;
north=0;
speed_intervals=[0,5,10,15,20,25,30,35];
%speed_intervals=[4,25];
%speed_intervals=[0,20,30];
nSpeedIntervals=length(speed_intervals);

[indices_dir,lim_dir]=sector_bins(T.(var_dir_ref),nSectors,north);
[indices_speed,lim_speed] = speed_bins(T.(var_speed_ref),speed_intervals);
 
% m_number=zeros(nSectors+1,nSpeedIntervals);
% m_mean=zeros(nSectors+1,nSpeedIntervals);
% m_mode=zeros(nSectors+1,nSpeedIntervals);
% m_mspeed=zeros(nSectors+1,nSpeedIntervals);
% m_TI_kelly=zeros(nSectors+1,nSpeedIntervals);
%  


cd 'Shear_Plots'
for s=1:nSectors+1
    for sp=1:nSpeedIntervals

        indices=intersect(indices_dir{s},indices_speed{sp});
     
        mU=zeros(1,length([z1:1:z2]));
        stdU=zeros(1,length([z1:1:z2]));
        lenU=zeros(1,length([z1:1:z2]));
        hU=heights;
     
%         for ih=1:length(mU)/2
%             var1=strcat(boom_dir1,'_WS_',num2str(heights(ih)),'m_Avg');
%             var2=strcat(boom_dir2,'_WS_',num2str(heights(ih)),'m_Avg');
%             ave=(T.(var1)+T.(var2))/2;
%             mU(ih)=mean(ave(indices));
%             stdU(ih)=std(ave(indices));
%             lenU(ih)=length(indices);                           
%         end

        for ih=1:length(heights)
            var1=strcat(boom_dir,'_WS_',num2str(heights(ih)),'m_Avg');
            var2=strcat(boom_dir,'_WS_',num2str(heights(ih)),'m_Avg');
            ave=(T.(var1)+T.(var2))/2;
            tt(1:length(indices),ih)=ave(indices);
        end
    
       %distribution of alpha 
        alpha_list=[];
        alpha_list_weight=[];
 
        for i=1:length(indices)
            fit=polyfit(log(heights([z1:z2])),log(tt(i,z1:z2)),1);
            weight=T.weight_months(i)*T.weight_days(i);
            if abs_on
                fit(1)=abs(fit(1));
            end
            alpha_list=[alpha_list;fit(1)];
            alpha_list_weight=[alpha_list_weight,transpose(weight)];
        end
       
        alpha_list=real(alpha_list);
        alpha_list_weight=transpose(alpha_list_weight);
    
        if length(alpha_list)<size_crit
            alpha_mean=nan;
            alpha_stderror=nan;
            alpha_50=nan;
            alpha_90=nan;
            alpha_mode=nan;
            spmean=mean(tt(1:length(indices),zref));
        else
           
            iinan=~isnan(alpha_list);
            alpha_list=alpha_list(iinan);
            alpha_list_weight=alpha_list_weight(iinan);
            
            dist_alpha=alpha_list;
%             alpha_mean=nanmean(alpha_list);
%             alpha_stderror=nanstd(alpha_list)/sqrt(length(alpha_list));

            alpha_mean=mean(alpha_list.*alpha_list_weight);
            alpha_stderror=std(alpha_list.*alpha_list_weight)/sqrt(length(alpha_list));

            spmean=nanmean(tt(1:length(indices),zref));
            
      
            %alpha_50=interp1(linspace(0.5/length(alpha_list), 1-0.5/length(alpha_list), length(alpha_list))', sort(alpha_list), 50*0.01, 'spline');
            %alpha_90=interp1(linspace(0.5/length(alpha_list), 1-0.5/length(alpha_list), length(alpha_list))', sort(alpha_list), 90*0.01, 'spline');
       
            out=wprctile(alpha_list,[50,90],alpha_list_weight,6);
            alpha_50=out(1);
            alpha_90=out(2);
            
            %plotting
            fig=figure;
            [pdfw,~, intervals] = pdf_weighted(alpha_list,alpha_list_weight,[-0.5:0.0125:0.5]);
            hb=bar(intervals, pdfw);
            hb.FaceColor = [0.3010 0.7450 0.9330];
            [maxi,id]=max(pdfw);
            mode=intervals(id);
            
            %fig=figure;
           %hist=histogram(alpha_list,[-0.5:0.0125:0.5],'Normalization','pdf');
%             
            

            %ind_max=find(max(hist.Values)==hist.Values);
            %mode=max((hist.BinEdges(ind_max)+hist.BinEdges(ind_max+1))/2);
            
            limy=maxi;
            
            alpha_mode=mode;
             hold on
             plot(alpha_mean*[1 1],[0,limy],'b')
             plot(alpha_50*[1 1],[0,limy],'r')
             plot(alpha_90*[1 1],[0,limy],'g') 
             plot(alpha_mode*[1 1],[0,limy],'k')  
             xlim([-0.5,0.5])
             %ylim([0,12])
             leg_data=strcat('data ','{ }',num2str(heights(z1)),'m-',num2str(heights(z2)),'m');
             leg_alpha=strcat('\alpha_{ave}=',num2str(round(alpha_mean,2)));
             leg_alpha_50=strcat('\alpha_{50}=',num2str(round(alpha_50,2)));
             leg_alpha_90=strcat('\alpha_{90}=',(num2str(round(alpha_90,2))));
             leg_alpha_mode=strcat('\alpha_{mode}=',(num2str(round(alpha_mode,2))));
      
             xlabel('shear exponent \alpha [-]','FontSize',fs)
             ylabel('pdf [-]','FontSize',fs)
             
            
             title(strcat('Met Mast Wind Shear in Sector:',{' '},char(lim_dir(s)),{' '},'and speed interval:',{' '},char(lim_speed(sp)) ));
             legend(leg_data,leg_alpha,leg_alpha_50,leg_alpha_90,leg_alpha_mode,'FontSize',fs,'Location','northwest')
             name=strcat(print_label,'_Met Mast Wind Shear_',num2str(heights(z1)),'m-',num2str(heights(z2)),'m_',num2str(s),'oo',num2str(nSectors),'_',num2str(sp),'oo',num2str(nSpeedIntervals),'.png');
             saveas(fig,name)
             hold off 
       end
       
      m_number(s,sp)=length(alpha_list);
      m_mean_r(s,sp)=round(alpha_mean,2);
      m_stderror_r(s,sp)=round(alpha_stderror,4);
      m_mode_r(s,sp)=round(alpha_mode,2);
      m_mspeed(s,sp)=spmean;
      m_mean(s,sp)=alpha_mean;
      m_mode(s,sp)=alpha_mode;
      
      m_TI_kelly(s,sp)=alpha_mode/(1+4*(alpha_mean-alpha_mode));
    end
end

cd '..'

%tt=table(transpose(lim_dir),m_mean_r(:,1),m_mean_r(:,2),m_mean_r(:,3),m_mean_r(:,4),'VariableNames',{'Sector',lim_speed{:}})

tt=array2table(m_mean_r,'VariableNames',{lim_speed{:}});
tt.('Sector')=transpose(lim_dir);
tt = [tt(:,end) tt(:,1:end-1)]
writetable(tt,strcat(print_label,'_MetMast_shear_mean_',num2str(heights(z1)),'m_',num2str(heights(z2)),'m.xlsx'),'WriteRowNames',true) 

tt2=array2table(m_stderror_r,'VariableNames',{lim_speed{:}})
tt2.('Sector')=transpose(lim_dir)
tt2 = [tt2(:,end) tt2(:,1:end-1)]
writetable(tt2,strcat(print_label,'_MetMast_shear_stderror_',num2str(heights(z1)),'m_',num2str(heights(z2)),'m.xlsx'),'WriteRowNames',true) 


%tt=table(transpose(lim_dir),m_TI_kelly(:,1),m_TI_kelly(:,2),m_TI_kelly(:,3),m_TI_kelly(:,4),'VariableNames',{'Sector',lim_speed{1},lim_speed{2},lim_speed{3},lim_speed{4}})
%writetable(tt,strcat(print_label,'_MetMast_shear_TI_Kelly_',num2str(heights(z1)),'m_',num2str(heights(z2)),'m.xlsx'),'WriteRowNames',true)


function [histw, vinterval] = histwc(vv, ww, nbins)
  minV  = min(vv);
  maxV  = max(vv);
  delta = (maxV-minV)/nbins;
  vinterval = linspace(minV, maxV, nbins)-delta/2.0;
  histw = zeros(nbins, 1);
  for i=1:length(vv)
    ind = find(vinterval < vv(i), 1, 'last' );
    if ~isempty(ind)
      histw(ind) = histw(ind) + ww(i);
    end
  end
end

function [histw, vinterval] = histwcv(vv, ww, nbins)  
  minV  = min(vv);
  maxV  = max(vv);
  delta = (maxV-minV)/nbins;
  vinterval = linspace(minV, maxV, nbins)-delta/2.0;
  histw = zeros(nbins, 1);
  indX  = arrayfun(@(xx) find(vinterval < vv(xx), 1, 'last'), 1:length(vv));
  arrayfun(@(xx) evalin('caller', ['histw(indX(', sprintf('%d', xx),')) = histw(indX(', sprintf('%d', xx),')) + ww(', sprintf('%d', xx),');']), 1:length(vv));
end
