function m_HT(H,T)
% modified by jwo to to overlay correctly
% Produce HT-scatter diagrams and THmax (and H-MWD diagram if directional)

% For dev:
% Considerations on directional extremes (HFH)
% Consider appropriatenes of zero-crossing
        
% Filename
filename = strcat(H.name,'_HT_',H.label,'_',T.label,'_',H.ttt_str); filename(isspace(filename)) = '_';
filename2 = strcat(H.name,'_THmax_',H.label,'_',T.label,'_',H.ttt_str); filename2(isspace(filename2)) = '_';
if H.subseries
    if strcmp(H.group{2},'Directional')
        subheading = strcat({' '},H.dirstruc.label,{' ['},H.dirstruc.unit,{']: '});
    else
        H.sublabel = '';
        H.subunit ='-';
        subheading = '';
    end
else
    subheading = '';
end


% Check time step and duration (1 second precision)
t_step = unique(round(24*3600*(H.time(2:end)-H.time(1:end-1))))/3600; % in hours
if length(t_step)  > 1, warning('Time series is non-equidistant - Check Nyears!'),end

% Initialize
N = length(H.time); NaNs = 0;
HT_Scatter = zeros(length(T.bins)-1,length(H.bins)-1,length(H.title));
Hmax = zeros(N,length(H.title)); THmax = zeros(N,length(H.title));
Hsum = zeros(length(H.bins)-1,1); Tsum = zeros(length(H.bins)-1,1);

% For each timestep
if isfield(H,'subtype')
    workspace_name  = strcat(H.name,'_',H.legend,'_',H.ttt_str,'_HT_workspace_',H.subtype,'.mat'); %  load('HT_workspace_directional.mat');
else
    workspace_name  = strcat(H.name,'_',H.legend,'_',H.ttt_str,'_HT_workspace.mat'); %  load('HT_workspace.mat');
end
 
if ~exist(workspace_name,'file')
    
h = waitbar(0,['HT: ' datestr(H.time(1),'yyyy-mm-dd HH:MM:SS')]);

for i = 1:N
    
    waitbar(i/N,h,['HT: ' datestr(H.time(i),'yyyy-mm-dd HH:MM:SS')]);
    % Gather HT - directional bins
    if H.subseries
        try
            dirbin = H.dirbin(i);
        catch
            dirbin = H.mthbin(i); % PDG: dirty fix: use dirbin for mthbin...
        end
    else
        dirbin = 1; 
    end
    
    if ~isempty(H.data(i).H)
        % Find index of H's and match with T's
        [n,xout]=histc(H.data(i).H,H.bins);
        HTtmp = [T.data(i).T,xout];
        % Make HT scatter table
        HTtable = zeros(length(T.bins)-1,length(H.bins)-1);
        for j = 1:length(H.bins)-1
            tmp = histc(HTtmp(find(HTtmp(:,2)==j),1),T.bins);
            % Add wave to H bin j (without last bin from histc)
            HTtable(:,j) = HTtable(:,j) + reshape(tmp(1:length(T.bins)-1),[],1);
            % Hsum and Tsum for each H.bin
            Hsum(j) = Hsum(j) + sum(H.data(i).H(find(H.data(i).H >= H.bins(j) & H.data(i).H < H.bins(j+1))));
            Tsum(j) = Tsum(j) + sum(T.data(i).T(find(H.data(i).H >= H.bins(j) & H.data(i).H < H.bins(j+1))));
        end
        HT_Scatter(:,:,dirbin) = HT_Scatter(:,:,dirbin) + HTtable;
        
        % Save Hmax and associated T in each sea state
        Hmax(i,dirbin) = max(H.data(i).H);
        rH = find(Hmax(i,dirbin) == H.data(i).H);
        THmax(i,dirbin) = T.data(i).T(rH);
        
        % Omni max
        Hmax(i,length(H.title)) = max(H.data(i).H);
        rH = find(Hmax(i,length(H.title)) == H.data(i).H);
        THmax(i,length(H.title)) = T.data(i).T(rH);
        
    else
        NaNs = NaNs+1; % Count time steps with NaNs
    end
    
end
close(h)

HT_Scatter(:,:,length(H.title)) = sum(HT_Scatter(:,:,:),3);
Hsum(end+1) = sum(Hsum);
Hcumsum = cumsum(Hsum(1:end-1))./Hsum(end);
Tsum(end+1) = sum(Tsum);
% Shift omni to first !!!
HT_Scatter = circshift(HT_Scatter,[0 0 1]);
Hmax = circshift(Hmax,[0 1]);
THmax = circshift(THmax,[0 1]);

else
    load(workspace_name)
end

% Sum omni into new matrix
% Convert to percent
% nHT = sum(sum(HT_Scatter(:,:,1)));
% P = HT_Scatter(:,:,:)/nHT*100;
% Convert to waves per year
% nYears = (H.time(end)-H.time(1))/365.25;
nYears = (N-NaNs)*H.ttt(4)/60/24/365.25; % Do not count NaNs!!!
P = HT_Scatter(:,:,:)/nYears;
% Add totals:
P(end+1,:,:) = sum(P(:,:,:),1);
P(:,end+1,:) = sum(P(:,:,:),2);
% Add Cumulated:
P(end+1,:,:) = cumsum(P(end,:,:));
P(:,end+1,:) = cumsum(P(:,end,:));
P(end,end-1:end,:) = 0; P(end-1:end,end:end,:) = 0;
% Create binlabels for scatter tables
for l = 1:length(H.bins)-1,  x_str(l) = {['[',num2str(H.bins(l)),'-',num2str(H.bins(l+1)),'[']}; end
for l = 1:length(T.bins)-1,  y_str(l) = {['[',num2str(T.bins(l)),'-',num2str(T.bins(l+1)),'[']}; end

if isfield(H,'subtype')
    save(workspace_name); %  load('HT_workspace_directional.mat');
else
    save(workspace_name); %  load('HT_workspace.mat');
end

% Print and Tables
for i = 1:length(H.title)
    
     % Contour - save('Tyne_Tees_35years','HT_Scatter(:,:,i),'nYears')
%     Hlab = H.bins(2:end) - (H.bins(2:end)-H.bins(1:end-1))/2;
%     Tlab = T.bins(2:end) - (T.bins(2:end)-T.bins(1:end-1))/2;
%     levels = [10^6 3*10^5 10^5 3*10^4 10^4 3*10^3 10^3 3*10^2 10^2 3*10^1 10^1];
%     [C,h] = contour(Tlab,Hlab,HT_Scatter(:,:,i)'/nYears,levels); hold on, grid on; clabel(C,h);
%     
%     set(gca, 'ylim', [0 14]), set(gca, 'xlim', [0 25]), ylabel('H(m)'), xlabel('T(s)');
%     title(['Tyne Tees',{'35years'}]); legend(['HT, N = ' num2str(round(sum(sum(HT_Scatter(:,:,i)'/nYears)))) ])
%     print(gcf,'-dpng','-r600','Tyne_Tees_35years.png'), close
    
    % HT scatter (without inf)
    HTS = log10(HT_Scatter(:,:,i)); HTS(HTS==0) = 0.1; HTS(isinf(HTS)) = 0;
    
    % Contour
    Hlab = H.bins(2:end) - (H.bins(2:end)-H.bins(1:end-1))/2;
    Tlab = T.bins(2:end) - (T.bins(2:end)-T.bins(1:end-1))/2;
    [C,h]=contour(Hlab,Tlab,HTS,10),hold on, grid on; colormap(H.ColorMap)
    legend = {['Contours of ln(N), N = ' num2str(P(end-1,end-1,1),'%2.1f')]};
    %     ntimes = 4; % NB: function 'interp2' requires toolbox (for smoothing)
    %     Hlab = interp(H.bins(2:end)-(H.bins(2:end)-H.bins(1:end-1))/2,2^ntimes); Hlab = Hlab(1:end-2^ntimes+1);
    %     Tlab = interp(T.bins(2:end)-(T.bins(2:end)-T.bins(1:end-1))/2,2^ntimes); Tlab = Tlab(1:end-2^ntimes+1);
    %     contour(Hlab,Tlab,interp2(HTS,ntimes),10), hold on, grid on, legend('Smoothed contours of ln(N_{HT})')
    set(gca, 'xlim', [H.bins(1) H.bins(end)]), set(gca, 'ylim', [T.bins(1) T.bins(end)])
    xlabel({[H.label ' (' H.unit ')']}), ylabel({[T.label ' (' T.unit ')']});
    set(gca,'ytick',T.bins), set(gca,'xtick',1:1:max(H.bins)),  set(gca,'xtickLabel',1:1:max(H.bins)), m_xticklabel_rotate([],45,[]);
    title([H.name,strcat('HT',H.ttt_str,subheading,{' '},H.title{i})],'FontWeight','normal');
    %         Xlab = H.bins(2:end) - (H.bins(2:end)-H.bins(1:end-1))/2;
    %         Ylab = T.bins(2:end) - (T.bins(2:end)-T.bins(1:end-1))/2;
    %         contour(Xlab,Ylab,PS,10);
    
    % Add Goda steepness (2010 formulation)
% jwo
% A = 0.18;
%     slope = 1/1000;
%     beta = -1.5*pi*abs(H.xyz(3))*(1+11*(tan(slope*pi/180))^(4/3));
%     for m = 1:length(T.bins)
%         L0(m) = 9.82/(2*pi)*T.bins(m)^2;
%         Hb(m) = L0(m)*A*(1-exp(beta/L0(m)));
%     end
%     plot(Hb,T.bins,'color',H.ColorOrder(4,:),'LineStyle','-')
%     legend_Goda = {['Breaking limit cf. Goda (A = ' num2str(A,'%2.2f') ')']};
    k0 = find(Hcumsum<0.99,1,'last');
    % Find mean in each Hmax bin!
    for k = 1:length(H.bins)-1
        Iavg = find(Hmax(:,i) >= H.bins(k) & Hmax(:,i) < H.bins(k+1));
        if k<k0
            Havg(k) = nan;
            Tavg(k) = nan;
        elseif ~isempty(Iavg)
            Havg(k) = nanmean(Hmax(Iavg,i));  Hstd(k) = nanstd(Hmax(Iavg,i));
            Tavg(k) = nanmean(THmax(Iavg,i)); Tstd(k) = nanstd(THmax(Iavg,i));
        else
            Havg(k) = nan;
            Tavg(k) = nan;
        end
        clear Iavg
    end
    % Remove NaNs
    Havg = Havg(isfinite(Havg)); Hstd = Hstd(isfinite(Hstd));
    Tavg = Tavg(isfinite(Tavg)); Tstd = Tstd(isfinite(Tstd));
    
    %plot(Havg,Tavg,'color',H.ColorOrder(2,:),'linestyle','none','marker','o')
    %legend_avg = [{['Mean +/- 2 x std. dev. of T_{Hmax}']}];
    
    %[estimates,fval,exitflag,model] = m_fit(Havg,Tavg,'Poly'); [sse, FittedCurve] = model(estimates);
    % jwo
    % plot(Havg,FittedCurve,'color',H.ColorOrder(2,:))
    %legend_fit = {['Least-Square fit: T_{Hmax} = ' num2str(estimates(1),'%2.2f') 'x' 'H_{max}' strcat('^{',num2str(estimates(2),'%2.2f'),'}')]};
    
    % Plot dots
    for j = 1:length(H.bins)-1, for k = 1:length(T.bins)-1, if HTS(k,j) > 0
                plot((H.bins(j+1)+H.bins(j))/2,(T.bins(k+1)+T.bins(k))/2,'ko','MarkerSize',max(HTS(k,j),1)), hold on, end, end, end
    legend_dots = {'Weighted occurrence'};
    
% jwo
% Plot +/-2xstd.dev. - round two - to get legend right!
%     for k = 1:length(Havg)
%         plot([     Havg(k)      Havg(k)],[Tavg(k)+2*Tstd(k) Tavg(k)-2*Tstd(k)],'color',H.ColorOrder(2,:))
%         plot([0.99*Havg(k) 1.01*Havg(k)],[Tavg(k)+2*Tstd(k) Tavg(k)+2*Tstd(k)],'color',H.ColorOrder(2,:))
%         plot([0.99*Havg(k) 1.01*Havg(k)],[Tavg(k)-2*Tstd(k) Tavg(k)-2*Tstd(k)],'color',H.ColorOrder(2,:))
%     end
%     
%     % Size and Position
%     legend([legend_contour legend_Goda legend_avg legend_fit legend_dots])
    fpos = get(gcf,'Position'); set(gcf,'Position',[fpos(1)   fpos(2)   fpos(3)   fpos(4)],'PaperPositionMode','auto')
    set(gca,'Position',[0.12 0.14 0.73 0.76])
    ylabel_pos = get(get(gca,'ylabel'),'Position'); set(get(gca,'ylabel'),'position',[-0.09             ylabel_pos(2) ylabel_pos(3)]);
    xlabel_pos = get(get(gca,'xlabel'),'Position'); set(get(gca,'xlabel'),'position',[xlabel_pos(1) 0.7*xlabel_pos(2) xlabel_pos(3)]);
    print(gcf,'-dpng',H.reso,strcat(filename,'_',num2str(i-1,'%02.0f'),'_',H.title{i},H.figure)), close
%     m_table(strcat(filename,'_',num2str(i-1,'%02.0f'),'_',H.title{i}),...
%         [H.name,strcat('Average annual number of waves ',H.ttt_str,subheading,{' - '},H.title{i})],...
%         {[H.label ' (' H.unit ')']},{[T.label ' (' T.unit ')']},...
%         [x_str 'Total' 'Accum'],[y_str 'Total' 'Accum'],P(:,:,i)','Precision',1,'Width',55,H.ascii,H.table,'ColorMap',H.ColorMap);
    m_table(strcat(filename,'_',num2str(i-1,'%02.0f'),'_',H.title{i}),...
        [H.name,strcat('Average annual number of waves ',H.ttt_str,subheading,H.title{i})],...
        {[H.label ' (' H.unit ')']},{[T.label ' (' T.unit ')']},...
        [x_str 'Total' 'Accum'],[y_str 'Total' 'Accum'],P(:,:,i)','Precision',1,'Width',55,H.ascii,H.table,'ColorMap',H.ColorMap);
    
    
    %% Plot THmax vs Hmax
%     if i == 1 % Omni only
        subplot(2,2,1:2),hold on, grid on
        legendstr = 'Max. H in seastate and associated T';
        plot(Hmax(:,i),THmax(:,i),'.k','MarkerSize',2),hold on
        
        % Find parameters in conditional distribution of Tp
%         [row,col,H2] = find(Hmax(:,i));  H2 = reshape(H2,1,[]); % only non-zero values for fit
%         [row,col,T2] = find(THmax(:,i)); T2 = reshape(T2,1,[]); % only non-zero values for fit
        IHT = find(Hmax(:,i)>=0.01*H.bins(end-1) & THmax(:,i)>=0); % Allow cut-off with respect to minimum Hmax for fit
        H2 = reshape( Hmax(IHT,i),1,[]); % only non-zero values for fit
        T2 = reshape(THmax(IHT,i),1,[]); % only non-zero values for fit
        hs = H.bins(2:end) - (H.bins(2:end)-H.bins(1:end-1))/2;
%         muhat = NaN(1,length(hs)); ssqrhat = NaN(1,length(hs));
        for j = find(min(H2)>H.bins,1,'last'):length(H.bins)-1
            ttmp = T2(find(H2 >= H.bins(j) & H2 < H.bins(j+1)));
            n = length(ttmp);
            muhat(j) = sum(log(ttmp))/n;
            ssqrhat(j) = sum((log(ttmp)-muhat(j)).^2)/n;
        end
        ii = find(~isnan(muhat));
        [a_est(i,:),fval,exitflag,model] = m_fit(hs(ii),muhat(ii),'FitMean');  [sse, aCurve] = model(a_est(i,:));
        [b_est(i,:),fval,exitflag,model] = m_fit(hs(ii),ssqrhat(ii),'FitVar'); [sse, bCurve] = model(b_est(i,:));
%         Tcurve.mean  = a_est(i,:); % output added FLD (Tcurve = m_HT_fld(H,T))
%         Tcurve.stedv = b_est(i,:); % output added FLD PDG: Should be as varargout or alternative through the structures (H and/or T)

        % Plot mean
        subplot(2,2,3),hold on, grid on
        plot(hs(ii),muhat(ii),'+b'),hold on,grid on
        plot(hs(ii),aCurve,'r')
        legend('Observed \mu^{ }',['\mu^{} = ',...
            num2str(a_est(i,1),'%5.4f'),' + ',...
            num2str(a_est(i,2),'%5.4f'),' \times ','H_{max}','^{',...
            num2str(a_est(i,3),'%5.4f'),'}'],'Location','NorthOutside') % PDG: 4 decimals needed for accuracy!
        xlabel('H_{max}  [m]'),ylabel('\mu'),xlim([H.bins(1) H.bins(end)])
        set(gca,'xtick',H.bins), set(gca,'xtickLabel',H.bins), m_xticklabel_rotate([],45,[]);
        
        % Plot var
        subplot(2,2,4),hold on, grid on
        plot(hs(ii),ssqrhat(ii),'+b'),hold on,grid on
        plot(hs(ii),bCurve,'r')
        legend('Observed \sigma^2',['\sigma^2 = ',...
            num2str(b_est(i,1),'%5.4f'),' + ',...
            num2str(b_est(i,2),'%5.4f'),' \times exp(-',...
            num2str(b_est(i,3),'%5.4f'),'','H_{max}',')'],'Location','NorthOutside') % PDG: 4 decimals needed for accuracy!
        xlabel('H_{max}  [m]'),ylabel('\sigma^2'),xlim([H.bins(1) H.bins(end)])
        set(gca,'xtick',H.bins), set(gca,'xtickLabel',H.bins), m_xticklabel_rotate([],45,[]);
        
        % Plot fraktions
        farve = [H.ColorOrder(4,:);H.ColorOrder(3,:);H.ColorOrder(2,:);H.ColorOrder(3,:);H.ColorOrder(4,:)];
        ltype = cellstr(strvcat('-.','--','-','--','-.')); P2=[0.05,0.25,0.5,0.75,0.95];
        Tp = nan(length(H.bins)-1,length(P2));
        
        for j = 1:length(P2)
            Tp((1:length(aCurve))+find(min(H2)>H.bins,1,'last')-1,j)   = invlognorm(P2(j),aCurve,bCurve);
            subplot(2,2,1:2), plot(hs(ii),Tp((1:length(aCurve))+find(min(H2)>H.bins,1,'last')-1,j),char(ltype(j)),'Color',farve(j,:),'LineWidth',1)
            legendstr = strvcat(legendstr,['T_{Hmax,',int2str(P2(j)*100),'%}']);
        end
        
        % Plot DNV
        subplot(2,2,1:2), plot(H.bins,2.94*H.bins.^0.5,'color',H.ColorOrder(5,:),'LineWidth',1.0), grid on
        legendstr = strvcat(legendstr,'T_{Hmax} = 2.94\times H_{max}^{0.5} (DNV RP-C205)');
        
        % Annotate
        set(gca,'ytick',T.bins), set(gca,'xtick',H.bins), set(gca,'xtickLabel',H.bins), m_xticklabel_rotate([],45,[]);
        legend(legendstr,'Location','SouthEastOutside'), xlabel('H_{max} [m]'); ylabel(['T_{Hmax} [s]']),xlim([H.bins(1) H.bins(end)]),ylim([T.bins(1) T.bins(end)])
        xlabel_pos = get(get(gca,'xlabel'),'Position'); set(get(gca,'xlabel'),'position',[xlabel_pos(1) 0.7*xlabel_pos(2) xlabel_pos(3)]);
        set(gca,'Position',[0.13 0.49 0.41 0.40])
        title([H.name,strcat('Associated T',H.ttt_str,{' '},H.title{i})],'FontWeight','normal');
        
        % Print
        print(gcf,'-dpng',H.reso,strcat(filename2,'_',num2str(i-1,'%02.0f'),'_',H.title{i},H.figure)), close
%         for k = 1:length(P2), xbin_str{k} = num2str(P2(k)); end
%         for k = 1:(length(H.bins)-1), ybin_str{k} = num2str(hs(k)); end
%         m_table(strcat(filename2,'_',num2str(i-1,'%02.0f'),'_',H.title{i}),...
%             [H.name,strcat('THmax ',H.ttt_str,{' - '},H.title{i})],...
%             {'Quantile'},{strcat(H.label,' (',H.unit,')')},...
%             xbin_str,ybin_str(ii),Tp','Precision',4,'Width',55,H.ascii,H.table);

%         end

end

if isfield(H,'subtype')
    filename2 = strcat(filename2,'_',H.subtype);
end
 m_table(strcat(filename2,'_params'),...
            [H.name,strcat('THmax ',H.ttt_str,{' - params'})],...
            {'param'},{'value'},...
            [{'a1'} {'a2'} {'a3'} {'b1'} {'b2'} {'b3'}],H.title,[a_est b_est]','Precision',4,'Width',55,H.ascii,H.table);

% % Calc of THmax for known Hmax
% % lognorm.inv in Excel takes std.dev (=variance^0.5) as input, whereas invlognorm in matlab takes the variance
% P = [0.05 0.25 0.5 0.75 0.95];
% Hmax = [9.1 11 11.7 12.7 13.5 14.3];
% mu=a(1)+a(2)*Hmax.^a(3);
% ssqr = b(1) + b(2)*exp(-b(3)*Hmax);
% for i=1:length(P)
%     THmax(i,:) = invlognorm(P(i),mu,ssqr);
% end
% disp(num2str(THmax,'%5.1f\t'))


%% MWD table
if H.subseries
    % Init (with omni i first column)
    HDIR = zeros(size(P,2)-1,size(P,3));
    for i = 1:size(P,3)
        HDIR(:,i) = sum(P(1:end-2,1:end-1,i),1)';
    end
    HDIR = circshift(HDIR,[0 size(P,3)-1]); % Shift omni to last
    HDIR = [Hsum./(HDIR(1:end,end)*nYears) HDIR]; % Add Havg column at start
    HDIR(:,end+1) = [Tsum./(HDIR(1:end,end)*nYears)]; % Add Tavg column at end
    % Table
    filename3 = strcat(H.name,'_HT_',H.label,'_',H.dirstruc.label,'_',H.ttt_str); filename3(isspace(filename3)) = '_';
    m_table(filename3,[{H.name} {['Average annual number of waves ',H.ttt_str]}],...
        {[H.dirstruc.label ' (' H.dirstruc.unit ')']},{[H.label ' (' H.unit ')']},...
        ['H_{avg}' H.title(2:end) H.title(1) 'T_{avg}'],[x_str 'Total'],HDIR',...
        'Precision',1,'Width',55,H.ascii,H.table); % ,'ColorMap',H.ColorMap
end

end

