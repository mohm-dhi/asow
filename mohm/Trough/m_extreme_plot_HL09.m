function [EVA] = m_extreme_plot(EVA)
%
% Plots results of unconstrained and constrained directional distribution fits
%
% Revision History
% HFH 2011-10-27:	First release
% HFH 2011-11-07:   Added truncated Weibull dist (ifit=2)
% HFH 2011-12-06:   Added ndec and xtick fields for user control of XTick
% HFH 2012-07-30:   Linear lnN(mode) relationship by default
% HFH 2012-08-20:   Added output of maximum distribution (for use with m_optimize_Tr)
% HFH 2012-09-26:   Added dims field allowing directions to be skipped from
% HFH 2012-10-05:   Corrected calculation of return periods in annual max. series
% HFH 2013-04-04:   Added proper plotting of GPD function (ifit=6)
% HFH 2014-05-09:   Added handling of Glukhovskiy distribution
% PDG 2014-07-16:   Added Hessian confidence limits (Gumbel only) (HFH input)
% HFH 2014-08-15:   Corrected error in convolution when EVA.sqrfac~=1
% PDG 2014-08-26:   Modified output to full EVA structure
% PDG 2014-09-14:   Corrected error in Hessian confidence limits (relevant for Tr>5yr mainly) (HFH input)
% BJE 2016-11-14:   Major reworking of the input to the m_table
% functionality, corrected for confidence limit and seasonal/directional
% output
% HFH 2018-10-08:   Added automatic subplot generation when more than 16 subsectors and changed
% plotflag to assume vector input (see m_extreme.m header for effect of different flags)
% XABR 2019-09-24:  Number of digits to display increased to 4 for the
% fitting parameters (location,shape,scale)
%
% PDG: Filename
heading = strcat({'Extreme '},EVA.label,{' '},EVA.legend,EVA.ttt_str);
filename = strcat(EVA.name,'_Extreme_',EVA.item,'_',EVA.legend,EVA.ttt_str,'_',...
    EVA.EVtype,'_',num2str(EVA.EVcrit(end),'%2.2f'),'_',...
    num2str(EVA.EVdist),'_SF=',num2str(EVA.sqrfac),'_',EVA.estimationmethod); % PDG added: ,'_SF=',num2str(EVA.sqrfac) %     'IET=',num2str(EVA.intereventtime,'%2.1f'),'h_','IEL=',num2str(EVA.intereventlevel,'%2.1f'),'_',...    % SJA removed as filename was too long

if isfield(EVA,'ShortTermDist'), filename = strcat(filename,'_',EVA.ShortTermDist); end
if strcmp(EVA.item(1),'C') && ~strcmp(EVA.unit,'m/s'), filename = strcat(filename,'_',EVA.vref); end % PDG for Crest
filename(isspace(filename)) = '_';
if length(EVA.group) > 1
    if strcmp(EVA.group(2),'Directional')
        if isfield(EVA,'ShortTermDist')
            filename = strcat(filename,'_Directional_(',EVA.itemInfo.Hm0.dirstruc.label,')'); heading = strcat(heading,' Directional (',EVA.itemInfo.Hm0.dirstruc.label,')');
        else
            filename = strcat(filename,'_Directional_(',EVA.sublabel,')'); heading = strcat(heading,' Directional (',EVA.sublabel,')');
        end
    end
    if strcmp(EVA.group(2),'Monthly'),     filename = strcat(filename,'_Monthly');     heading = strcat(heading,' Monthly');end
    if strcmp(EVA.group(2),'Seasonal'),    filename = strcat(filename,'_Seasonal');    heading = strcat(heading,' Seasonal');end
    if strcmp(EVA.group(2),'Monsoon'),     filename = strcat(filename,'_Monsoon');     heading = strcat(heading,' Monsoon');end
    if strcmp(EVA.group(2),'Flood_Ebb'),   filename = strcat(filename,'_Flood_Ebb');   heading = strcat(heading,' Flood_Ebb');end
end
if isfield(EVA,'optstr'), filename = strcat(filename,EVA.optstr); end % JOAB
%
% only use Gumbel if EVtype=='AMP';
if strcmpi(EVA.EVtype,'AMP') && EVA.EVdist~=4; warning('Only Gumbel (ifit==4) should be used together with AMP selection method'); end
%
% determine number of dimensions
EVA.ndim = size(EVA.X,2)-1;
%
if isfield(EVA,'ShortTermDist') % These are mode data. Prepare for convolution
    if ~isfield(EVA,'lnN')
        % Calculate an appropriate value of lnN (used in convolution)
        tmp = sortrows(EVA.mode_N_shape(:,1:2),-1);
        
        % Fit linear relationsship to the modes and ln(avg(N)-STD(N)) data
        % from the historical storms
        if size(tmp,1)>100
            % bin data in 10 bins (10 largest modes and the remaining data in
            % equally spaced bins)
            stp = [linspace(tmp(end,1),tmp(10,1),10),2*tmp(1,1)];
        else
            % bin data in equally spaced bins of approx. 10 pts/bin
            stp = linspace(tmp(end,1),tmp(1,1)+1e-10,ceil(size(tmp,1)/10)+1);% stp(end+1) = 2*tmp(1,1);
        end
        [~,out]=histc(tmp(:,1),stp);
        mode_lnN = NaN(length(stp)-1,2);
        for i=1:length(stp)-1;
            mode_lnN(i,:) = [mean(tmp(  out==i,1)),log(mean(tmp(out==i,2))-std(tmp(out==i,2)))]; % PDG: this fails (imaginary) if std and/or mean is large!
        end
        % lnN_params=polyfit(mode_lnN(~isnan(mode_lnN(:,1)),1),mode_lnN(~isnan(mode_lnN(:,1)),2),1); % estimate parameters of linear relationship
        % PDG: Replacing above by added account for imaginary values and for 1 fitting point only (0 orcer hotfix!)
        Ireal = [];for k = 1:length(mode_lnN), if isreal(mode_lnN(k,:)) && ~isnan(sum(mode_lnN(k,:))), Ireal = [Ireal k]; end, end
        if length(Ireal) > 1
            lnN_params=polyfit(mode_lnN(Ireal,1),mode_lnN(Ireal,2),1); % estimate parameters of linear relationship
        else
            lnN_params=polyfit(mode_lnN(Ireal,1),mode_lnN(Ireal,2),0); % estimate parameters of linear relationship
            lnN_params(2) = lnN_params; lnN_params(1) = 0;
        end

        lnN_format='linear'; % flag defining format in which lnN information comes
        if ismember(4,EVA.plotflag)
            close all,plot(tmp(:,1),log(tmp(:,2)),'.'),hold on,grid on
            plot(mode_lnN(:,1),mode_lnN(:,2),'+r')
            plot([min(tmp(:,1)),max(tmp(:,1))],...
                lnN_params(1)*[min(tmp(:,1)),max(tmp(:,1))]+lnN_params(2),'r')
            xlabel('Mode [m]'),ylabel('ln(N)')
            print('-dpng','-r300',strcat(filename,'_lnN_f.png')), close
        end
        
    else
        lnN        = EVA.lnN;
        lnN_format = 'user-defined';
    end
    % PDG: This should be taken/adopted from m_Short_Term_Dist to avoid double
    shape_format = 'constant';
    switch EVA.ShortTermDist
        case 'Forristall_H'
            shape = 2.126;
        case 'Forristall_C'
            if ~isfield(EVA,'shape')
                % we also compute the shape parameter from data
                %(however, no subtraction of 1 STD)
                tmp   = sortrows(EVA.mode_N_shape(:,[1,3]),-1);
                tmp   = tmp(1:round(0.25*size(tmp,1)),2);
                shape = mean(tmp); clear tmp
            else
                shape = EVA.shape;
            end
        case 'Rayleigh_H' % PDG: added 
            shape = 2.000;
        case 'Rayleigh_T'  % FLD: added 
            shape = 2.000;
        case 'Naess_H' % PDG: added
            shape = 2.000;
        case 'Glukhovskiy_H'
            if ~isfield(EVA,'shape')
                % we also compute the shape parameter from data
                tmp          = sortrows(EVA.mode_N_shape(:,[1,3]),-1);
                shape_params = polyfit(tmp(:,1),tmp(:,2),1);
                shape_format = 'linear';
                if ismember(4,EVA.plotflag)
                    close all,plot(tmp(:,1),tmp(:,2),'.'),hold on,grid on
                    plot([min(tmp(:,1)),max(tmp(:,1))],...
                        shape_params(1)*[min(tmp(:,1)),max(tmp(:,1))]+shape_params(2),'r')
                    xlabel('Mode [m]'),ylabel('Glukhovskiy shape [-]')
                    try
                        print('-dpng','-r300',strcat(filename,'_G_shape.png')), close
                    end
                end
                
            else
                shape = EVA.shape;
            end
        case 'Battjes_Groenendijk_H'
            % discontinous shape with discontinuity at Htr
            Htr = (0.35+5.8*tan(EVA.slope))*EVA.WaterDepth;
            shape_format = 'B&G2000';
        case 'Weibull_H'
            shape = EVA.alpha;
            
    end
end
%
% misc. plot control
% % Make default title strings for all subplots if EVA.divisions is undefined
% if isempty(EVA.divisions)
%     for i=1:EVA.ndim; EVA.divisions{i} = ['Sector ',int2str(i)]; end;
%     EVA.divisions=[EVA.divisions,'Omni'];
% end


if EVA.ndim==0
    sb = [1,1,1];
    fs0 = 8;
elseif EVA.ndim<=2
    sb = [3,1,3];
    fs0 = 6;
elseif EVA.ndim<=4
    sb = [4,2,5];
    fs0 = 6;
elseif EVA.ndim<=8
    sb = [6,2,9];
    fs0 = 5.2;
elseif EVA.ndim<=12
    sb = [6,3,13];
    fs0 = 4;
elseif EVA.ndim<=16
    sb = [6,4,17];
    fs0 = 4;
else
    sb = [ceil(EVA.ndim/4)+1,4,ceil(EVA.ndim/4)*4+1];
    fs0 = 15/ceil(EVA.ndim/4); % make font size dependent on number of subplots
end
fs0 = 8; fs = fs0; % PDG

% Only plot the directional plots without the omni-direction JOAB
if ismember(6,EVA.plotflag) && size(EVA.dims,2)>2 
    EVA.dims = EVA.dims(1:end-1);
    EVA.ndim = length(EVA.dims);
    sb(3) = sb(3)-1;
    switch sb(1)
        case 3
            sb(1) = sb(1)-1;
        case 4
            sb(1) = sb(1)-1; 
        otherwise
            sb(1) = sb(1)-1; 
    end
end

% if strcmp(EVA.label,'H') || strcmp(EVA.label,'C'), EVA.label = strcat(EVA.label,'_{mp}'); end
ylbl = 'T_R [years]'; % PDG
%


% % define x-tik mark location for all plots
% xmx   = ceil(1.1*max(max(EVA.X(:,end)))^EVA.sqrfac);
% dx = round(xmx/8*10)/10;
% if isempty(EVA.xtick)
%     xtick = 0:dx:dx*8;
%     %         if EVA.sqrfac<1; xtick=[xtick(1),xtick(3:end)]; end % PDG: ?
% else
%     xtick = EVA.xtick;
% end
% textpos = floor(min(min(EVA.Location)))+xtick(1)^(1/EVA.sqrfac) + [...
%     0.05*(xtick(end)^(1/EVA.sqrfac)-floor(min(min(EVA.Location))))+xtick(1)^(1/EVA.sqrfac),...
%     linspace(0.35*(xtick(end)^(1/EVA.sqrfac)-floor(min(min(EVA.Location))))+xtick(1)^(1/EVA.sqrfac),0.9*(xtick(end)^(1/EVA.sqrfac)-floor(min(min(EVA.Location))))+xtick(1)^(1/EVA.sqrfac),length(EVA.T))]; % PDG
%
% % y-tick marks
% ytik = [1:1:10,20:10:100,200:100:1000,2e3:1e3:1e4,2e4:1e4:1e5];
% i0   = find(ytik<EVA.T(1),1,'last'); if isempty(i0); i0=1; end
% i1   = find(ytik>EVA.T(end),1,'first');
% ytik = ytik(i0:i1);
% ytikmajor = 10.^(0:6); ytikmajor = ytikmajor(ytikmajor<=EVA.T(end));
% ymn = -4; %-log(-log(ytik(1)));%floor(-log(-log(0.1)));
% ymx = ceil(-log(-log(1-1/EVA.T(end))));
% %
% % define number of decimals for plots
% if isempty(EVA.ndec)
%     ndec  = [0 2;1 1;100 0]; % intervaller that defines ndec
%     ndec  = ndec(find(ndec(:,1)<=xmx,1,'last'),2); % define ndec based on value of xmx
% else
%     ndec = EVA.ndec;
% end
% fmstr = cellstr(['%',int2str(ndec+4),'.',int2str(ndec),'f']);
%
% set default confidence limits
if isempty(EVA.ConfLimits); EVA.ConfLimits  = [0.025,0.975]; end
%

% Hessian confidence intervals (normz = 0.5 gives central estimate, normal dist., i: return period, j: directions)
if EVA.Hessian && EVA.EVdist == 4 % (Gumbel only)
    if isempty(EVA.ConstrainedParams)
        mu    = [EVA.UnconstrainedParams(1,end,1)];
        sigma = [EVA.UnconstrainedParams(2,end,1)];
    else
        mu    = [EVA.ConstrainedParams(1,:,1)]; % ,EVA.UnconstrainedParams(1,end,1)
        sigma = [EVA.ConstrainedParams(2,:,1)]; % ,EVA.UnconstrainedParams(2,end,1)
    end
    try
        EVA.R_LB = ReturnGumbelConfInterval(EVA.T,EVA.X,mu,sigma,-1.96).^EVA.sqrfac; % 2.5%
        EVA.R_UB = ReturnGumbelConfInterval(EVA.T,EVA.X,mu,sigma, 1.96).^EVA.sqrfac; % 97.5%
    catch
        EVA.R_LB = NaN(size(EVA.T,2),size(EVA.X,2));
        EVA.R_UB = NaN(size(EVA.T,2),size(EVA.X,2));
        warning('Hessian conf intervals not found!') % PDG quick fix (34843) - in case of inappropriate Uncons params, typically as a result of bad constraining
    end
end

if EVA.EVdist==6
%     X_est = linspace(0,-EVA.UnconstrainedParams(2,end)/EVA.UnconstrainedParams(1,end)-eps,2000)'; 
    X_est = linspace(0,EVA.UnconstrainedParams(2,end)/abs(EVA.UnconstrainedParams(1,end))-eps,2000)'; % PDG: replaced above - relevant for positive shape parameter!
    X_cum = X_est;
elseif EVA.EVdist==7
    X_est = linspace(0,EVA.UnconstrainedParams(3,end)-eps,2000)';
    X_cum = X_est;
    %
    % if EVA.EVdist==5
    %     X_cum = linspace(min(EVA.GAM),3*max(max(EVA.X)),2000)';
    %     i_cum = find(X_cum>max(EVA.GAM));
else
    X_est = linspace(0,log10(EVA.T(end))*max(max(EVA.X(:,:,1))),2000)';
    X_cum = linspace(1e-3,log10(EVA.T(end))*max(max(EVA.X(:,:,1))),2000)'+EVA.Location(1,end,1);
end
%
% determine if constrained set of distributions exists
if ~isempty(EVA.ConstrainedParams)
    nfit=2;
else
    nfit=1;
end
P_cum = ones(length(X_cum),nfit);
%
% confidence limits
if ~isempty(EVA.ConfLimits) && EVA.N_BootStrap>0
    cflim = min(max(round(EVA.ConfLimits*EVA.N_BootStrap),1),EVA.N_BootStrap);
    if min(cflim)<10; disp(['WARNING m_extreme_plot: ',num2str(100*min(EVA.ConfLimits),'%4.1f'),...
            '% confidence limit not well established from only ',int2str(EVA.N_BootStrap),' bootstraps']) ; end
    if EVA.N_BootStrap-max(cflim)<10; disp(['WARNING m_extreme_plot: ',num2str(100*max(EVA.ConfLimits),'%4.1f'),...
            '% confidence limit not well established from only ',int2str(EVA.N_BootStrap),' bootstraps']) ; end
else
    cflim=[];
end
%
% loop through distribution sets and directions
for k=1:nfit
    
    switch k
        case 1
            Params = EVA.UnconstrainedParams;
            %            imx = EVA.ndim+1;
        case 2
            Params = EVA.ConstrainedParams;
            %            imx = EVA.ndim;
    end
    
    for i=EVA.dims  %EVA.ndim+1
        
        if EVA.EVdist==1
            error('ifit==1 is not valid')
            P_obs = 1-linspace(lambda*Nyrs-0.44,1-0.44,lambda*Nyrs)/(lambda*Nyrs+0.12); %% !!! Gringorten plotting position
            P0 = exp(-(g/b)^a);
            X_obs = EVA.X(:,i);
            i_est = find(X_cum>EVA.GAM(i));
            P_Xmp = 1 - exp(-((X_est-EVA.GAM(i))/b)^a);
            Rmp   = b*(-log(P0./(T*EVA.XLAM(i)))).^(1/a) + EVA.GAM(i);
            
        elseif EVA.EVdist==2  % Truncated Weibull
            
            if k==1
                try
                    %                     XP(1,i).obs = [sort(EVA.X(~isnan(EVA.X(:,i)),i)),...
                    %                         ((1:EVA.lambda(i)*EVA.Nyrs)' + 0.5/sqrt(Params(1,i))-0.6)/...
                    %                         (EVA.lambda(i)*EVA.Nyrs + 0.2 + 0.23/sqrt(Params(1,i)))];
                    XP(1,i).obs = [sort(EVA.X(~isnan(EVA.X(:,i)),i)),...
                        ((1:sum(~isnan(EVA.X(:,i))))' + 0.5/sqrt(Params(1,i))-0.6)/...
                        (sum(~isnan(EVA.X(:,i))) + 0.2 + 0.23/sqrt(Params(1,i)))];
                catch
                    disp('stop')
                    
                end
            end
            P0 = exp(-(EVA.Location(i)/Params(2,i)).^Params(1,i));
            XP(k,i).est  = [X_est,1 - 1/P0*exp(-(X_est/Params(2,i)).^Params(1,i))];
            XP(k,i).est(X_est<EVA.Location(i),:)=NaN;
            XP(k,i).cum  = [X_cum,1 - 1/P0*exp(-(X_cum/Params(2,i)).^Params(1,i))];
            XP(k,i).R_mp = Params(2,i)*(-log(P0./(EVA.lambda(i)*EVA.T))).^(1/Params(1,i));
            %OneYrValue(k,i) = Params(2,i)*(-log(P0./EVA.lambda(i)))^(1/Params(1,i));
            for j=2:EVA.N_BootStrap+1 % compute bootstrap RP values
                P0 = exp(-(EVA.Location(1,i,j)/Params(2,i,j)).^Params(1,i,j));
                XP(k,i).BootStrap(:,j-1) = Params(2,i,j) * ...
                    (-log(P0./(EVA.T*EVA.lambda(1,i,j)))).^(1/Params(1,i,j));
            end;
            
        elseif EVA.EVdist==3  % 2-p Weibull to excess
            
            if k==1
                try
                    %                     XP(1,i).obs = [sort(EVA.X(~isnan(EVA.X(:,i)),i)),...
                    %                         ((1:EVA.lambda(i)*EVA.Nyrs)' + 0.5/sqrt(Params(1,i))-0.6)/...
                    %                         (EVA.lambda(i)*EVA.Nyrs + 0.2 + 0.23/sqrt(Params(1,i)))];
                    XP(1,i).obs = [sort(EVA.X(~isnan(EVA.X(:,i)),i)),...
                        ((1:sum(~isnan(EVA.X(:,i))))' + 0.5/sqrt(Params(1,i))-0.6)/...
                        (sum(~isnan(EVA.X(:,i))) + 0.2 + 0.23/sqrt(Params(1,i)))];
                catch
                    disp('stop')
                    
                end
            end
            XP(k,i).est  = [X_est+EVA.Location(i),1 - exp(-(X_est/Params(2,i)).^Params(1,i))];
            XP(k,i).cum  = [X_cum,1 - exp(-((X_cum-EVA.Location(i))/Params(2,i)).^Params(1,i))];
            XP(k,i).R_mp = Params(2,i)*(-log(1./(EVA.lambda(i)*EVA.T))).^(1/Params(1,i)) + EVA.Location(i);
            if Params(1,i)>1
                XP(k,i).Mode = Params(2,i)*((Params(1,i)-1)/Params(1,i))^(1/Params(1,i));
            elseif Params(1,i)==1
                XP(k,i).Mode = 0;
            else
                XP(k,i).Mode = NaN;
            end
            %OneYrValue(k,i) = Params(2,i)*(-log(1./EVA.lambda(i)))^(1/Params(1,i)) + EVA.Location(i);
            for j=2:EVA.N_BootStrap+1 % compute bootstrap RP values
                XP(k,i).BootStrap(:,j-1) = Params(2,i,j)*(-log(1./...
                    (EVA.T*EVA.lambda(1,i,j)))).^(1/Params(1,i,j)) + EVA.Location(1,i,j);
            end;
        elseif EVA.EVdist==4 % Gumbel
            if k==1
                nonan = ~isnan(EVA.X(:,i)); % PDG
                XP(1,i).obs = [sort(EVA.X(nonan,i)),...
                    1-linspace(EVA.lambda(i)*EVA.Nyrs-0.44,1-0.44,...
                    round(EVA.lambda(i)*EVA.Nyrs))'/(EVA.lambda(i)*EVA.Nyrs+0.12)];
            end
            XP(k,i).est  = [X_est,exp(-exp(-(X_est-Params(1,i))/Params(2,i)))];
            XP(k,i).cum  = [X_cum,exp(-exp(-(X_cum-Params(1,i))/Params(2,i)))];
            XP(k,i).R_mp = Params(1,i) - Params(2,i)*log(-log(1-1./(EVA.lambda(i)*T2TA(EVA.EVtype,EVA.T))));
            XP(k,i).Mode = Params(1,i);
            for j=2:EVA.N_BootStrap+1 % compute bootstrap RP values
                XP(k,i).BootStrap(:,j-1) = Params(1,i,j) - ...
                    Params(2,i,j)*log(-log(1-1./(EVA.lambda(1,i,j)*T2TA(EVA.EVtype,EVA.T))));
                %XP(k,i).BootStrap(EVA.T==1,j-1) = Params(1,i,j);    % 1-year value is mode
            end
            
        elseif EVA.EVdist==5  % Exponential
            
            if k==1
                try
                    XP(1,i).obs = [sort(EVA.X(~isnan(EVA.X(:,i)),i)),...
                        ((1:sum(~isnan(EVA.X(:,i))))' + 0.5-0.6)/...
                        (sum(~isnan(EVA.X(:,i))) + 0.2 + 0.23)];
                catch
                    disp('stop')
                    
                end
            end
            XP(k,i).est  = [X_est+EVA.Location(i),1 - exp(-(X_est/Params(1,i)))];
            XP(k,i).cum  = [X_cum,1 - exp(-((X_cum-EVA.Location(i))/Params(1,i)))];
            XP(k,i).R_mp = Params(1,i)*(-log(1./(EVA.lambda(i)*EVA.T))) + EVA.Location(i);
            XP(k,i).Mode = 0;
            for j=2:EVA.N_BootStrap+1 % compute bootstrap RP values
                XP(k,i).BootStrap(:,j-1) = Params(1,i,j)*(-log(1./...
                    (EVA.T*EVA.lambda(1,i,j)))) + EVA.Location(1,i,j);
            end;
            
        elseif EVA.EVdist==6  % Generalized Pareto
            
            if k==1
                % Gringorten plotting position..
%                 XP(1,i).obs = [sort(EVA.X(:,i)),...
%                     1-linspace(EVA.lambda(i)*EVA.Nyrs-0.44,1-0.44,...
%                     EVA.lambda(i)*EVA.Nyrs)'/(EVA.lambda(i)*EVA.Nyrs+0.12)];
                XP(1,i).obs = [sort(EVA.X(:,i)),...
                    1-linspace(EVA.lambda(i)*EVA.Nyrs-0.44,1-0.44,...
                    round(EVA.lambda(i)*EVA.Nyrs))'/(EVA.lambda(i)*EVA.Nyrs+0.12)]; % PDG: added round as for Gumbel
            end
            XP(k,i).est  = [X_est+EVA.Location(i),1 - (1 + Params(1,i)/Params(2,i)*X_est).^(-1/Params(1,i))];
            XP(k,i).est(XP(k,i).est(:,1)-EVA.Location(i)>-Params(2,i)/Params(1,i),:) = NaN; % Values outside supported range
            XP(k,i).cum  = [X_cum,1 - (1 + Params(1,i)/Params(2,i)*(X_cum-EVA.Location(i))).^(-1/Params(1,i))];
            %XP(k,i).cum(XP(k,i).cum(:,1)-EVA.Location(i)>-Params(2,i)/Params(1,i),:) = NaN; % Values outside supported range
            XP(k,i).R_mp = Params(2,i)/Params(1,i)*((1./(EVA.lambda(i)*EVA.T)).^-Params(1,i) - 1) + EVA.Location(i);
            
            %x = m + s.*-log(1-1./(EVA.lambda(i)*EVA.T))
            
            %OneYrValue(k,i) = Params(2,i)/Params(1,i)*((1./EVA.lambda(i)).^-Params(1,i) - 1) + EVA.Location(i);
            XP(k,i).Mode = NaN; %!!!!!!!!!!
            
            for j=2:EVA.N_BootStrap+1 % compute bootstrap RP values
                XP(k,i).BootStrap(:,j-1) = Params(2,i,j)/Params(1,i,j) * ...
                    ((1./(EVA.lambda(1,i,j)*EVA.T)).^-Params(1,i,j) - 1) + EVA.Location(1,i,j);
            end;
            
            
        elseif EVA.EVdist==7  % Upper Bound Weibull
            
            if k==1
                XP(1,i).obs = [sort(EVA.X(~isnan(EVA.X(:,i)),i)),...
                    ((1:sum(~isnan(EVA.X(:,i))))' + 0.5/sqrt(Params(1,i))-0.6)/...
                    (sum(~isnan(EVA.X(:,i))) + 0.2 + 0.23/sqrt(Params(1,i)))];
            end
            UB0 = Params(3,end); % we use the omnidirectional upper bound for all directions
            
            %                         P_untrc=exp(-((UB-X)/b).^a);
            %                         P0= exp(-((UB-g)/b)^a);
            %                         %P_trunc=(P_untrc-P0)/(1-P0);
            %                         n_untrc = (1./(y*Method(mthd).lambda)*(1-P0)).^-1;
            %                         R_plot  = UB - b*exp(log(-log(1-1./n_untrc))/a);
            %                         n_untrc = (1./(Method(mthd).T*Method(mthd).lambda)*(1-P0)).^-1;
            %                         R_Xmp = UB - b*exp(log(-log(1-1./n_untrc))/a);
            %                         n_storms = a/b/(1-P0)*(((UB-X)/b).^(a-1)).*(exp(-((UB-X)/b).^a)) .* Xstp;
            %
            %                         P_untrc=exp(-((UB-Xcum)/b).^a);
            %                         P_trunc=(P_untrc-P0)/(1-P0);
            %                         ncum = [ncum,(1-P_trunc)*Method(mthd).lambda];
            
            P0 = exp(-((UB0-EVA.Location(i))/Params(2,i)).^Params(1,i));
            XP(k,i).est  = [X_est,(exp(-((UB0-X_est)/Params(2,i)).^Params(1,i))-P0)/(1-P0)];
            XP(k,i).est(X_est<EVA.Location(i),:)=NaN;
            XP(k,i).cum  = [X_cum,(exp(-((UB0-X_cum)/Params(2,i)).^Params(1,i))-P0)/(1-P0)];
            XP(k,i).R_mp = UB0 - Params(2,i)*(-log(1-1./(EVA.lambda(i)*EVA.T*(1-P0)))).^(1/Params(1,i));
            for j=2:EVA.N_BootStrap+1 % compute bootstrap RP values
                UB = Params(3,end,i);
                P0 = exp(-((UB-EVA.Location(1,i,j))/Params(2,i,j)).^Params(1,i,j));
                XP(k,i).BootStrap(:,j-1) = UB - ...
                    Params(2,i,j)*(-log(1-1./(EVA.lambda(1,i,j)*EVA.T*(1-P0)))).^(1/Params(1,i,j));
            end;
            
            XP(k,i).Mode = NaN; %!!!!!!!!!!
            
            clear UB P0
        end
        
        % set 1-year value to mode of distribution in the case of AMP series
        if strcmpi(EVA.EVtype,'AMP')
            XP(k,i).R_mp(EVA.T==1) = XP(k,i).Mode;
        end
        
        if isfield(XP,'BootStrap')
            XP(k,i).BootStrap = sort(XP(k,i).BootStrap,2); % sort ascending
        end
        if i<=EVA.ndim
            P_cum(:,k) = P_cum(:,k).*XP(k,i).cum(:,2); %.^EVA.XLAM(i));
        end
        
        %%%% convolution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if k==nfit && isfield(EVA,'ShortTermDist')
            if EVA.EVdist==6
%                 tmp = linspace(0,-EVA.UnconstrainedParams(2,i,1)/EVA.UnconstrainedParams(1,i,1)-eps,1e4)+EVA.Location(1,i,1);
                tmp = linspace(0,EVA.UnconstrainedParams(2,i,1)/abs(EVA.UnconstrainedParams(1,i,1))-eps,1e4)+EVA.Location(1,i,1); % PDG: replaced above - relevant for positive shape parameter!
                xstp = tmp(2:end)-tmp(1:end-1);
                x    = tmp(1:end-1)+(tmp(2)-tmp(1))/2;
            elseif EVA.EVdist==7
                tmp  = linspace(EVA.Location(1,i,1),UB0,1e4);
                xstp = tmp(2:end)-tmp(1:end-1);
                x    = tmp(1:end-1)+(tmp(2)-tmp(1))/2;
            else
                tmp  = logspace(0,log10(4*max(EVA.X(:,i,1))-EVA.Location(1,i,1)),1e4)-1;
                xstp = tmp(2:end)-tmp(1:end-1);
                x    = (cumsum(xstp)-xstp/2)+EVA.Location(1,i,1);
            end
            xmx  = linspace(EVA.Location(1,i,1),5*max(EVA.X(:,i,1)),1000)'; % PDG: added: '0.8*' on 'EVA.Location(1,i,1)' but doesnt solve the underlying issue
            
            
            %                         X0 = UB + 1 - logspace(0,log10(UB+1),100000); X0(end) = 0; X0 = sort(X0);
            %                         X0 = X0(find(X0<=UB & X0>g0));
            %                         Xstp = NaN(size(X0));
            %                         Xstp(2:end-1) = (X0(3:end) - X0(1:end-2))/2;
            %                         Xstp(1) = Xstp(2)/2; Xstp(end) = Xstp(end-1)/2;
            
            if strcmpi(lnN_format,'linear')
                % ln(N) as a function of mode estimated by linear
                % regression. Calculate ln(N) and expand to appropriate
                % array size
                lnN = meshgrid(lnN_params(1)*x.^EVA.sqrfac+lnN_params(2),1:length(xmx));
            end
            if strcmpi(shape_format,'linear')
                % shape as a function of mode estimated by linear
                % regression. Calculate shape and expand to appropriate
                % array size
                shape = meshgrid(shape_params(1)*x.^EVA.sqrfac+shape_params(2),1:length(xmx));
            end
            if strcmpi(shape_format,'B&G2000')
                % compute shape with discontinuity at Htr
                shape = meshgrid([repmat(2,1,sum(x.^EVA.sqrfac<Htr)),repmat(3.6,1,sum(x.^EVA.sqrfac>=Htr))],1:length(xmx));
            end
            
            n_storms = pdf_funcs(EVA.EVdist,x,EVA.Location(1,i,1),Params(:,i,1)).*xstp;
            % some code that corrects small numerical errors or causes the
            % script to stop if the error is to large
            fac = 1/sum(n_storms); if abs(1-fac)>1e-2; warning(['Something wrong with n_storms, pls check!!! Direction ',int2str(i)]);end
            %if 1-sum(n_storms)<0; n_storms = n_storms * fac; end
            n_storms(1) = 1-sum(n_storms(2:end))-eps; % fix first value in n_storms so integral of pdf == 1
            
            P_xmax = exp(-exp(-lnN.*(((xmx)*((x).^-1)).^(shape*EVA.sqrfac)-1))) * n_storms'; %.^EVA.sqrfac !?! 
            XP(k,i).max = [xmx,P_xmax];
            
            % we do the same for the bootstrap distributions (if any)
            if EVA.N_BootStrap>0
                XP(k,i).BootStrapMax = NaN(length(xmx),EVA.N_BootStrap);    % init array
                cc = 0;
                for j=2:EVA.N_BootStrap+1                                   % compute bootstrap RP values
                    try
                        n_storms = pdf_funcs(EVA.EVdist,x,EVA.Location(1,i,j),Params(:,i,j)).*xstp;
                        fac = 1/sum(n_storms); if abs(1-fac)>1e-3; error('noget galt med n_storms');end
                        n_storms(1) = 1-sum(n_storms(2:end))-eps; % fix first value in n_storms so integral of pdf == 1
                        P_xmax = exp(-exp(-lnN.*(((xmx)*((x).^-1)).^shape-1))) * n_storms';
                        XP(k,i).BootStrapMax(:,j-1) = P_xmax;
                    catch
                        cc=cc+1;
                        warning(['Convolution of boostrap no ',int2str(j-1),' for dir. ',int2str(i),...
                            ' failed... Now failed ',int2str(cc),' of ',int2str(EVA.N_BootStrap),' for this direction.'])
                    end
                end
                XP(k,i).BootStrapMax = sort(XP(k,i).BootStrapMax,2,'descend'); % sort descending
            end
        end
    end
    
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R_mp = NaN(length(EVA.T),EVA.ndim+2,nfit,length(cflim)+1); % Initialize array for storing return period values
R_mp(:,1,:,:) = ndgrid(EVA.T,1,1:nfit,1:length(cflim)+1);      % Put return periods in first column
% R_mp(return periods, [directional dimensions or seasons with last==omni], [1==unconstraint; 2==
% constraint], central estimate & confidence limits)
R_mx = R_mp;
for i=EVA.dims  %EVA.ndim+1
    for k=1:nfit
        R_mp(:,i+1,k,1) = XP(k,i).R_mp.^EVA.sqrfac;
        if k==nfit && isfield(EVA,'ShortTermDist')
            R_mx(:,i+1,k,1) = interp1q(XP(k,i).max(:,2),XP(k,i).max(:,1),1-1./(EVA.lambda(i)*T2TA(EVA.EVtype,reshape(EVA.T,[],1)))).^EVA.sqrfac;
        end
        for jj=1:length(cflim)
            R_mp(:,i+1,k,jj+1) = XP(k,i).BootStrap(:,cflim(jj)).^EVA.sqrfac;
            if k==nfit && isfield(EVA,'ShortTermDist')
                R_mx(:,i+1,k,jj+1) = interp1q(XP(k,i).BootStrapMax(:,2),XP(k,i).max(:,1),1-1./(EVA.lambda(i)*T2TA(EVA.EVtype,reshape(EVA.T,[],1)))).^EVA.sqrfac;
            end
        end
    end
end

            
% MGO: get rid of imaginary parts in the plot. PDG: moved from below for loops to reduce repetition and get into outputs...
if ~isreal(R_mp)
    R_mp = real(R_mp); % BJE modified sqrt((real(R_mp(r,i+1,nfit,1))^2+imag(R_mp(r,i+1,nfit,1))^2)); % JOAB modified to utilize "real"
    R_mx = real(R_mx); % JOAB added R_mx directional check.
    warning('Some or all return period values are complex. Check the EVA output!'); % BJE added
end

% PDG: Added to output:
EVA.R_mp = R_mp;
EVA.R_mx = R_mx;
EVA.XP = XP;

if ismember(0,EVA.plotflag), return, end % HFH, 2018-10-08: Now, adding 0 to plotflag vector will overrule all other plotflags and completely abort plotting
clf
legend_str = [{['Data Point (N = ' num2str(size(EVA.X,1)) ')']},{'Unconstrained Fit'}]; % PDG
%
% define number of decimals for plots
if isempty(EVA.ndec)
    ndec = max(0,ceil(2-log10((max(max(EVA.X(:,:,1)))^EVA.sqrfac))));
    EVA.ndec = ndec; EVA.ndec; % PDG
else
    ndec = EVA.ndec;
end
ndec2 = ndec-1; % ndec for X-tickmarklabels

% fmstr1 = cellstr(['%',int2str(ndec+4),'.',int2str(ndec),'f']);
% fmstr2 = cellstr(	);
%
% % define x-tik mark location for all plots
if isempty(EVA.xtick)
    if isfield(EVA,'ShortTermDist')
        xmx = ceil(1.05*max(max(R_mx(:,2:end,nfit,1)))*10^(ndec-0))/10^(ndec-0); % PDG: added 5% and changed -2 to -0
    else
        xmx = ceil(1.05*max(max(R_mp(:,2:end,nfit,1)))*10^(ndec-0))/10^(ndec-0); % PDG: added 5% and changed -2 to -0
    end
    xmn = floor(min(min(EVA.X(:,:,1)))^EVA.sqrfac*10^(ndec-0))/10^(ndec-0); %and changed -2 to -0
    
    dx = [2,5,10]'*logspace(-5,5,11); dx = dx(:);
    [~,idx] = min(abs(8-(xmx-xmn)./dx));
    xtick = xmn:dx(idx):xmx; if xtick(end)<xmx; xtick(end+1)=xtick(end)+dx(idx); end
    clear xmn xmx dx idx
else
    xtick = EVA.xtick;
    % check that ndec is sufficiently large to not round of xtick
    cc=0;
    while sum(abs(round(xtick*10^ndec2)-xtick*10^ndec2))>1e-9
        ndec2=ndec2+1;
        if ndec<ndec2;
            ndec = ndec2;
            warning('Number of decimals (ndec) increased to %i to accomodate user-defined xtick\n',ndec)
        end
        cc=cc+1;
        if cc==5; error('Something wrong with user-defined xtick. Use more round values for xticks'); end
    end
end
textpos = xtick(1) + (xtick(end)-xtick(1))*[0.05,linspace(0.35,0.9,length(EVA.T))];
% textpos = [0.05*(xtick(end)^(1/EVA.sqrfac))+xtick(1)^(1/EVA.sqrfac),...
%     linspace(0.35*(xtick(end)^(1/EVA.sqrfac))+xtick(1)^(1/EVA.sqrfac),...
%     0.9*(xtick(end)^(1/EVA.sqrfac))+xtick(1)^(1/EVA.sqrfac),length(EVA.T))];
% textpos = floor(min(min(EVA.Location)))+xtick(1)^(1/EVA.sqrfac) + [...
%     0.05*(xtick(end)^(1/EVA.sqrfac)-floor(min(min(EVA.Location))))+xtick(1)^(1/EVA.sqrfac),...
%     linspace(0.35*(xtick(end)^(1/EVA.sqrfac)-floor(min(min(EVA.Location))))+xtick(1)^(1/EVA.sqrfac),0.9*(xtick(end)^(1/EVA.sqrfac)-floor(min(min(EVA.Location))))+xtick(1)^(1/EVA.sqrfac),length(EVA.T))]; % PDG

%
% y-tick marks
ytik = [1:1:10,20:10:100,200:100:1000,2e3:1e3:1e4,2e4:1e4:1e5,2e5:1e5:1e6,2e6:1e6:1e7];
i0   = find(ytik>=EVA.T(1)/EVA.lambda(end),1,'first'); if isempty(i0); i0=1; end
i1   = find(ytik>EVA.T(end),1,'first');
ytik = ytik(i0:i1);
ytikmajor = 10.^(0:6); ytikmajor = ytikmajor(ytikmajor<=EVA.T(end) & ytikmajor>=ytik(1));
% ymn = -4; %-log(-log(ytik(1)));%floor(-log(-log(0.1)));
% ymx = ceil(-log(-log(1-1/EVA.T(end))));
%ymn = XP(1,end).obs(1,2)-1; %(XP(1,end).obs(2,2)-XP(1,end).obs(1,2));
if isempty(EVA.ylim)
    ymx = -log(-log(1-1./(EVA.lambda(end)*ytik(end))));
    %     ymn = ymx-7; %!!!!!!!!!!!!
    ymn = -log(-log(XP(1,end).obs(1,2))); % PDG (not valid for multidim POT)
else
    ymn = EVA.ylim(1);
    ymx = EVA.ylim(2);
end

neg = 1; % PDG: invert XTickLabel and vals for low water level
if strcmpi(EVA.item,'WL_low') 
    neg = -1;
    R_mp(:,2) = neg*R_mp(:,2);
end
%
%
% sb=[4,3,12];
%

for i=EVA.dims
    if i<=EVA.ndim
        subplot(sb(1),sb(2),i),hold on, grid on
        if EVA.rep_dir
            tstr = EVA.title(EVA.rep_dirs(i)+1);
        else
            tstr = EVA.title(i+1);
        end
    else
        subplot(sb(1),sb(2),sb(3):sb(1)*sb(2)),hold on, grid on;
        tstr = EVA.title(1);
    end
    %         kmx = 2;
    %         cum = 0;
    %     else
    %         %         fs = 1.3*fs0;
    %
    %         kmx = 1;
    %         cum = 2;
    %     end
    
    % plot data points (peaks)
    plot(XP(1,i).obs(:,1).^EVA.sqrfac,-log(-log(XP(1,i).obs(:,2))),'color',EVA.ColorOrder(1,:),'LineStyle','none','Marker','o','MarkerSize',3,'MarkerFaceColor',EVA.ColorOrder(1,:))
    
    % plot cumulative probability distributions - unconstrained
    plot(XP(1,i).est(:,1).^EVA.sqrfac,-log(-log(XP(1,i).est(:,2))),'color',EVA.ColorOrder(1,:),'LineStyle','-');
    
    % plot cumulative pcrobability distributions - constrained
    if nfit==2 %EVA.ConstrainFit && EVA.ndim > 1
        plot(XP(2,i).est(:,1).^EVA.sqrfac,-log(-log(XP(2,i).est(:,2))),'color',EVA.ColorOrder(2,:),'LineStyle','-');
        if i==EVA.ndim+1; legend_str = [legend_str,strcat('Constrained Fit (f_{const.}= ',num2str(EVA.constfac,'%2.0f'),')')]; end
    end
    
    % plot accumulated probability (omni-plot only)
    if i==EVA.ndim+1 && length(EVA.title)>1 && ismember(5,EVA.plotflag) % HFH: add 5 to plotflag vector for plot of product of constrained dists
        %         for k=1:cum; plot(XP(1,i).cum(:,1),-log(-log(P_cum(:,k))),'color',EVA.ColorOrder(k,:),'LineStyle','-'); end
        plot(XP(1,i).cum(:,1).^EVA.sqrfac,-log(-log(P_cum(:,1))),'color',EVA.ColorOrder(4,:),'LineStyle','-');
        legend_str = [legend_str,{'Accumulated Unconstrained Fits'}];
        if length(EVA.title)>1
            plot(XP(1,i).cum(:,1).^EVA.sqrfac,-log(-log(P_cum(:,2))),'color',EVA.ColorOrder(5,:),'LineStyle','-');
            legend_str = [legend_str,{'Accumulated Constrained Fit'}];
        end
    end
    
    % plot distribution of maximum wave/crest
    if isfield(EVA,'ShortTermDist')
        plot(XP(nfit,i).max(:,1).^EVA.sqrfac,-log(-log(XP(nfit,i).max(:,2))),'color',EVA.ColorOrder(3,:),'LineStyle','-');
        if i==EVA.ndim+1; legend_str = [legend_str,{['Convolution w. short-term dist. (' EVA.label(1) '_{max})']}]; end
    end
    
    % plot confidence limits
    if ~isempty(cflim)
        for j = cflim%([1 length(cflim)]) % PDG
            if isfield(XP(k,i),'BootStrapMax'); % plot confidence limits of maximum wave/crest
%                 plot(XP(k,i).max(:,1).^EVA.sqrfac,-log(-log(XP(k,i).BootStrapMax(:,j))),'color',EVA.ColorOrder(3,:),'LineStyle','--');
                if nfit==2
                    plot(XP(k,i).max(:,1).^EVA.sqrfac,-log(-log(XP(k,i).BootStrapMax(:,j))),'color',EVA.ColorOrder(2,:),'LineStyle','--');
                else
                    plot(XP(k,i).max(:,1).^EVA.sqrfac,-log(-log(XP(k,i).BootStrapMax(:,j))),'color',EVA.ColorOrder(1,:),'LineStyle','--');
                end
            else
                if EVA.ConstrainFit
                    plot(XP(nfit,i).BootStrap(:,j).^EVA.sqrfac,-log(-log(1-1./(EVA.lambda(i)*T2TA(EVA.EVtype,EVA.T)))),'color',EVA.ColorOrder(2,:),'LineStyle','--','LineWidth',1)
                else
                    plot(XP(nfit,i).BootStrap(:,j).^EVA.sqrfac,-log(-log(1-1./(EVA.lambda(i)*T2TA(EVA.EVtype,EVA.T)))),'color',EVA.ColorOrder(1,:),'LineStyle','--','LineWidth',1)
                end
            end
        end
        if i==EVA.ndim+1; legend_str = [legend_str,{'Bootstrap conf. lim. (Upper)'}]; end
        if i==EVA.ndim+1; legend_str = [legend_str,{'Bootstrap conf. lim. (Lower)'}]; end
    end
    
    % plot Hessian limits
    if EVA.Hessian && EVA.EVdist == 4
%         for j = cflim%([1 length(cflim)]) % PDG
%             if isfield(XP(k,i),'BootStrapMax'); % plot confidence limits of maximum wave/crest
%                 plot(XP(k,i).max(:,1).^EVA.sqrfac,-log(-log(XP(k,i).BootStrapMax(:,j))),'color',EVA.ColorOrder(3,:),'LineStyle','--');
%             else
%                 if EVA.ConstrainFit
                    plot(EVA.R_UB(:,i),-log(-log(1-1./(EVA.lambda(i)*T2TA(EVA.EVtype,EVA.T)))),'color',EVA.ColorOrder(3,:),'LineStyle','--','LineWidth',1)
                    plot(EVA.R_LB(:,i),-log(-log(1-1./(EVA.lambda(i)*T2TA(EVA.EVtype,EVA.T)))),'color',EVA.ColorOrder(3,:),'LineStyle','--','LineWidth',1)
%                 else
%                     plot(XP(nfit,i).BootStrap(:,j).^EVA.sqrfac,-log(-log(1-1./(EVA.lambda(i)*T2TA(EVA.EVtype,EVA.T)))),'color',EVA.ColorOrder(1,:),'LineStyle','--','LineWidth',1)
%                 end
%             end
%         end
        if i==EVA.ndim+1; legend_str = [legend_str,{'Hessian conf. lim. (Upper)'}]; end
        if i==EVA.ndim+1; legend_str = [legend_str,{'Hessian conf. lim. (Lower)'}]; end
    end
    
    %     % plot confidence limits of maximum wave/crest
    %     if ~isempty(cflim)
    %         for j = cflim([1 length(cflim)]) % PDG
    %             if isfield(XP(k,i),'BootStrapMax');
    %                 plot(XP(k,i).max(:,1),-log(-log(XP(k,i).BootStrapMax(:,j))),'color',EVA.ColorOrder(3,:),'LineStyle','--');
    %             end
    %         end
    %     end
    
    % print return period values on plots
    text(textpos(1),0.93*(ymx-ymn)+ymn,'T_R [years]','FontSize',fs) % ,'FontWeight','bold','BackGround','w'
    text(textpos(1),0.74*(ymx-ymn)+ymn,[EVA.label,'  [',EVA.unit,']'],'FontSize',fs) % ,'BackGround','w'
    if isfield(EVA,'ShortTermDist')
        text(textpos(1),0.53*(ymx-ymn)+ymn,strcat(EVA.label(1),'_{max}','  (',EVA.unit,')'),'FontSize',fs) % ,'BackGround','w'
    end
    
    for r = 1:length(EVA.T)
        text(textpos(r+1),0.88*(ymx-ymn)+ymn,int2str(EVA.T(r)),'FontSize',fs) % ,'BackGround','w','FontWeight','bold'        
        text(textpos(r+1),0.69*(ymx-ymn)+ymn,num2str(R_mp(r,i+1,nfit,1),['%',int2str(ndec+4),'.',int2str(ndec),'f']),'FontSize',fs) % ,'BackGround','w'
        if isfield(EVA,'ShortTermDist')
            text(textpos(r+1),0.53*(ymx-ymn)+ymn,num2str(R_mx(r,i+1,nfit,1),['%',int2str(ndec+4),'.',int2str(ndec),'f']),'FontSize',fs) % ,'BackGround','w'
        end
    end
    
     set(gca,'XTick',xtick,'XTickLabel',num2str(neg*xtick',['%',int2str(ndec+4),'.',int2str(max(0,ndec)),'f']),... % PDG: ndec2 -> ndec
        'XLim',xtick([1,end]),...
        'FontSize',fs,'YLim',[ymn ymx],...
        'YTick',unique(real(-log(-log(1-1./(EVA.lambda(i)*T2TA(EVA.EVtype,ytik)))))),'YTickLabel',' ')
    y = -log(-log(1-1./(EVA.lambda(i)*T2TA(EVA.EVtype,ytikmajor))));
    plot([xtick(1),xtick(end)],meshgrid(y,1:2),'color',EVA.ColorOrder(1,:),'LineStyle','--')
    text(xtick(1)*ones(size(ytikmajor))-0.02,y,int2str(ytikmajor'),'FontSize',fs,'HorizontalAlignment','right') % HFH correction 2012-10-04
    
    
    %     set(gca,'XTick',xtick.^(1/EVA.sqrfac),'XTickLabel',num2str(xtick',['%',int2str(ndec+4),'.',int2str(max(0,ndec-1)),'f']),...
    %         'XLim',xtick([1,end]).^(1/EVA.sqrfac),...
    %         'FontSize',fs,'YLim',[ymn ymx],...
    %         'YTick',-log(-log(1-1./(EVA.lambda(i)*T2TA(EVA.EVtype,ytik)))),'YTickLabel',' ')
    %     y = -log(-log(1-1./(EVA.lambda(i)*T2TA(EVA.EVtype,ytikmajor))));
    %     plot([xtick(1),xtick(end)].^(1/EVA.sqrfac),meshgrid(y,1:2),'color',EVA.ColorOrder(1,:),'LineStyle','--')
    %     text(floor(min(min(EVA.Location)))*ones(size(ytikmajor))-0.02,y,int2str(ytikmajor'),'FontSize',fs,'HorizontalAlignment','right') % HFH correction 2012-10-04
    %     if i==EVA.ndim+1
    %         xtick = get(gca,'XTick');
    %         xmx  = ceil(1.2*max(max(EVA.X(:,i)+EVA.Location(i)))^EVA.sqrfac);
    %         ndec = round(10/xmx);
    %         fmstr = cellstr(['%',int2str(ndec+4),'.',int2str(ndec),'f']);
    %         if EVA.sqrfac~=1
    %             xtick = round(tmp.^EVA.sqrfac*ndec)/ndec;
    %         end
    %     end
    title(tstr,'FontWeight','normal');
    % if kmx == 1; % PDG replacing line above
    %     title(EVA.title(1));
    % else
    %     title(EVA.title(i+1));
    % end
    ylabel(ylbl),
    xlabel([EVA.label,'  (',EVA.unit,')'])
    gpos = get(gca,'Position');
    if rem(i,sb(2)) == 0 && sb(3) == 3
        set(gca,'Position',[0.70*gpos(1) (1+(i-1)*0.05)*gpos(2) 1.1*gpos(3) 0.9*gpos(4)])
    elseif rem(i,sb(2)) == 0
        set(gca,'Position',[1.00*gpos(1) gpos(2) 1.2*gpos(3) 0.9*gpos(4)])
    elseif rem(i,sb(2)) == 1
        set(gca,'Position',[0.70*gpos(1) gpos(2) 1.2*gpos(3) 0.9*gpos(4)])
    elseif rem(i,sb(2)) == 2
        set(gca,'Position',[0.95*gpos(1) gpos(2) 1.2*gpos(3) 0.9*gpos(4)])
    end
    if i == EVA.ndim+1
        set(gca,'Position',[0.70*gpos(1) (1.5-sb(2)*0.11)*gpos(2) 1.1*gpos(3) 0.9*gpos(4)])
    end
end

% Special plotflag for MOOD!
if ismember(9,EVA.plotflag)
    title([strcat(EVA.name,EVA.xyz_str),heading],'FontWeight','normal');
    [legend_h, object_h, plot_h, text_strings] = legend([legend_str(1) 'Extreme Distribution' 'Confidence Limits'],'location','northeastoutside','FontSize',fs+1);
    a_pos = get(gca,'Position');
   
    % to replace numbers on top of plot:
%     A(1,1) = {'T_R [years]'};
%     A(1,2) = {[EVA.label,'  [',EVA.unit,']']};
%     for ii = 1:length(EVA.T)
%         A(ii+1,1) = {int2str(EVA.T(ii))};
%         A(ii+1,2) = {num2str(R_mp(ii,i+1,nfit,1),['%',int2str(ndec+4),'.',int2str(ndec),'f'])};
%     end
%             sprintf([blanks(1) A{1,1} blanks(5) A{1,2} '\n'... 
%                  blanks(8) A{2,1} blanks(16) A{2,2} '\n'...
%                  blanks(8) A{3,1} blanks(16) A{3,2} '\n'...
%                  blanks(7) A{4,1} blanks(16) A{4,2} '\n'...
%                  blanks(7) A{5,1} blanks(16) A{5,2} '\n'...
%                  blanks(6) A{6,1} blanks(16) A{6,2}]),...
    
    annotation('textbox',[legend_h.Position(1) a_pos(2) legend_h.Position(3) 0.28],'BackgroundColor','white','String',[...
        sprintf('Distribution:  Gumbel\n    Location = %2.2fm\n    Scale     = %2.2fm \n',EVA.UnconstrainedParams(1,end,1),EVA.UnconstrainedParams(2,end,1)),...
        sprintf('Events:  Annual Maxima\n'),...
        sprintf('Fitting:  Max. Likelihood\n'),...
        sprintf('Plotting pos.:  Gringorten\n'),...
        sprintf('Uncertainty:  Hessian\n\n'),...
        ],'FontSize',fs+1); % ,'FontName','FixedWidth', The ability to use the "\t" format in a listbox in MATLAB to place a tab in a string expression is not available
    
    print('-dpng',EVA.reso,strcat(filename,EVA.figure)), close
    
    % Add table
    for i = 1:length(EVA.T), ylabs(i) = {num2str(EVA.T(i))}; end
    table_row_title = {'Lower confidence','Central Estimate','Upper confidence'};
    m_table(filename,{strcat(filename,EVA.figure)},{'T_R [years]'},{[EVA.label,'  [',EVA.unit,']']},ylabs,table_row_title,...
        [EVA.R_LB R_mp(:,2:end) EVA.R_UB],EVA.ascii,EVA.table,'Precision',EVA.ndec);
    return
end

% Legend
[legend_h, object_h, plot_h, text_strings] = legend(legend_str,'location','northeastoutside');
% if EVA.ConstrainFit, set(plot_h(3),{'Color'},{'r'}); end
% if EVA.Hessian, set(plot_h(4),{'Color'},{'g'}); end

% Position
fpos = get(gcf,'Position');
set(gcf,'Position',[fpos(1)/6   fpos(2)/6   max(sb(2)*0.36,1.4)*fpos(3)   max(sb(1)*0.4,1.4)*fpos(4)],'PaperPositionMode','auto')
axes('position',[0,0,1,1],'visible','off');


% Title
text(0.5,0.96,[strcat(EVA.name,EVA.xyz_str),heading],'horizontalAlignment','center');

% Annotate
text(0.06,0.08,['Bias Correction: ',num2str(EVA.EVadju),'%'],'FontSize',fs);
text(0.06,0.06,['Event Selection: ',EVA.EVtype,' (',num2str(EVA.EVcrit(end),'%2.2f'),')'],'FontSize',fs);
text(0.06,0.04,['Inter-Event Time: ',num2str(EVA.intereventtime,'%2.1f'),'h'],'FontSize',fs);
text(0.06,0.02,['Inter-Event Level: ',num2str(EVA.intereventlevel,'%2.1f')],'FontSize',fs);

text(0.32,0.08,['Square factor = ' num2str(EVA.sqrfac)],'FontSize',fs);
text(0.32,0.06,['Dist. Type: ',EVA.LongTermDist],'FontSize',fs);
text(0.32,0.04,['Estimation Method: ',EVA.estimationmethod],'FontSize',fs);

% Bottom centered position of legend JOAB
if ismember(6,EVA.plotflag) && EVA.ndim>=2
    legend_h.Position = [0.42,0.18,0.17,0.05];
end

if isfield(EVA,'ShortTermDist')
    if strcmp(EVA.ShortTermDist,'Naess_H')
        text(0.32,0.02,['Short-Term Dist: ',EVA.ShortTermDist,' (\rho = ',num2str(EVA.rho),')'],'FontSize',fs);
    else
        text(0.32,0.02,['Short-Term Dist: ',EVA.ShortTermDist],'FontSize',fs);
    end
end

text(0.52,0.08,'Distribution Parameters: ','FontSize',fs);
if EVA.sqrfac~=1; unitstr = [' [',EVA.unit,']^{',num2str(1/EVA.sqrfac),'}']; else unitstr=EVA.unit; end
if ismember(EVA.EVdist,[2,3,6])
    text(0.52,0.06,'Location','FontSize',fs),text(0.59,0.06,['= ' num2str(EVA.Location(end,end,1),'%2.4f'),unitstr],'FontSize',fs);
    text(0.52,0.04,'Scale','FontSize',fs)   ,text(0.59,0.04,['= ',num2str(EVA.UnconstrainedParams(2,end,1),'%2.4e'),unitstr],'FontSize',fs);
    text(0.52,0.02,'Shape','FontSize',fs)   ,text(0.59,0.02,['= ',num2str(EVA.UnconstrainedParams(1,end,1),'%2.4f'),''],'FontSize',fs);
    if EVA.UnconstrainedParams(1,end,1) < 1 && EVA.EVdist~=6; warning('Omni-directional shape parameter is < 1'), end
elseif EVA.EVdist==4
    text(0.52,0.06,'Location','FontSize',fs),text(0.59,0.06,['= ' num2str(EVA.UnconstrainedParams(1,end,1),'%2.4f'),unitstr],'FontSize',fs);
    text(0.52,0.04,'Scale','FontSize',fs)   ,text(0.59,0.04,['= ',num2str(EVA.UnconstrainedParams(2,end,1),'%2.4e'),unitstr],'FontSize',fs);
elseif EVA.EVdist==5
    text(0.52,0.06,'Location','FontSize',fs),text(0.59,0.06,['= ' num2str(EVA.Location(end,end,1),'%2.4f'),unitstr],'FontSize',fs);
    text(0.52,0.04,'Scale','FontSize',fs)   ,text(0.59,0.04,['= ',num2str(EVA.UnconstrainedParams(1,end,1),'%2.4f'),unitstr],'FontSize',fs);
elseif EVA.EVdist==7
    text(0.52,0.06,'UB','FontSize',fs),text(0.59,0.06,['= ' num2str(UB0,'%2.2f'),unitstr],'FontSize',fs);
    text(0.52,0.04,'Scale','FontSize',fs)   ,text(0.59,0.04,['= ',num2str(EVA.UnconstrainedParams(2,end,1),'%2.4f'),unitstr],'FontSize',fs);
    text(0.52,0.02,'Shape','FontSize',fs)   ,text(0.59,0.02,['= ',num2str(EVA.UnconstrainedParams(1,end,1),'%2.4f'),''],'FontSize',fs);
    
end

if  EVA.ConstrainFit
    text(0.70,0.08,['Constraint Periods: ' num2str(EVA.Tcon) ' years'],'FontSize',fs);
end
if ~isempty(cflim);
    text(0.70,0.06,['Uncertainty Method: ' 'Bootstrap (',num2str(EVA.N_BootStrap) ')'],'FontSize',fs);
    text(0.70,0.04,['Confidence Limits (dashed lines): ' num2str(EVA.ConfLimits(1)*100) ' and ' num2str(EVA.ConfLimits(end)*100) '%'],'FontSize',fs);
end

filename = [filename '_CstFc' num2str(EVA.constfac,'%2.0f')];
print('-dpng',EVA.reso,strcat(filename,EVA.figure)), close

% % set x-axis
% for i=EVA.dims  %EVA.ndim+1
%     if i<=EVA.ndim; subplot(sb(1),sb(2),i)
%     else subplot(sb(1),sb(2),EVA.ndim+1:sb(1)*sb(2)); end
%     set(gca,'XTick',xtick.^(1/EVA.sqrfac),'XTickLabel',num2str(xtick',char(fmstr(1))),...
%         'XLim',[xtick(1),xtick(end)].^(1/EVA.sqrfac))
% end

% % Print/write table PDG: (% uncommented by fld)
subtitle{1} = 'All';
if ~isfield(EVA,'rep_dirs'), EVA.rep_dirs = 1:length(EVA.title); end
for i = 1:length(EVA.T), ylabs(i) = {num2str(EVA.T(i))}; end
if  isfield(EVA,'ShortTermDist') % length(XP) > 2 && isfield(XP(2,1),'max')
    % m_table(filename,{EVA.name subtitle{1}},{'MAX'},{'T_R'},circshift(EVA.title,[0 -1]),[ylabs ylabs],[R_mx(:,2:end,2)],EVA.ascii,EVA.table,'Precision',EVA.ndec); % R_mx(:,:,1) ;
    % m_table(filename,{EVA.name subtitle{1}},{[EVA.label(1) '_{max}' '  (' EVA.unit ')']},{'T_R [years]'},circshift(EVA.title,[0 -1]),ylabs,R_mx(:,2:end,2)',EVA.ascii,EVA.table,'Precision',EVA.ndec); % R_mx(:,:,1) ;
    try % constained values
        m_table(filename,{EVA.name subtitle{1}},{'T_R [years]'},{[EVA.label(1) '_{max}' '  (' EVA.unit ')']},ylabs,EVA.title([max(1,EVA.rep_dirs(end)-length(EVA.rep_dirs)) EVA.rep_dirs(1:end-1)+1]),circshift(R_mx(:,2:end,2),[0 1]),EVA.ascii,EVA.table,'Precision',EVA.ndec); % R_mx(:,:,1) ; % pdg rotated - and use 2 = constrained values!!!
    catch % else
        m_table(filename,{EVA.name subtitle{1}},{'T_R [years]'},{[EVA.label(1) '_{max}' '  (' EVA.unit ')']},ylabs,EVA.title([max(1,EVA.rep_dirs(end)-length(EVA.rep_dirs)) EVA.rep_dirs(1:end-1)+1]),circshift(R_mx(:,2:end,1),[0 1]),EVA.ascii,EVA.table,'Precision',EVA.ndec); % R_mx(:,:,1) ; % pdg rotated - and use 2 = constrained values!!!
    end
%     tmp = R_mx(:,2:end,1)'; save([filename '.mat'],'ylabs','tmp')
%     XP.R_mx = tmp;
    
% elseif length(XP) == 1 && isfield(XP(1,1),'max'); %added by fld
% %     m_table(filename,{EVA.name subtitle{1}},{'MAX'},{'T_R'},circshift(EVA.title,[0 -1]),[ylabs ylabs],[R_mx(:,2:end,2)],EVA.ascii,EVA.table,'Precision',EVA.ndec); % R_mx(:,:,1) ;
%     m_table(filename,{EVA.name subtitle{1}},{[EVA.label(1) '_{max}' '  (' EVA.unit ')']},{'T_R [years]'},circshift(EVA.title,[0 -1]),ylabs,R_mx(:,2:end,2)',EVA.ascii,EVA.table,'Precision',EVA.ndec); % R_mx(:,:,1) ;

else
    % BJE major reworking: this will now always write the constraint values
    % to table and all confidence limits as well directions etc
    % if omni only, unconstraint value is used...
    %
        if length(cflim) == 2
            % only two confidence limits + central estimate
            temp_legend = repmat({'Lower bnd','Central Est.','Upper bnd'},1,length(EVA.title)); % create a text string for the upper and lower bound values
            for m = 1:length(EVA.title)*3 table_row_title(m) = {[char(EVA.title(ceil(m/3))),...
                    ' ',char(temp_legend(m))]}; end % concat the confidence labels and the directional /seasonal labels to one 
            table_row_title = circshift(table_row_title,[0,(length(EVA.title)-1)*3]); %Make omni diectional label last rows by shifting the order of labels...
            cf_order = [2 1 3];
            
        elseif length(cflim) > 2
            % more than two (user defined) confidence limits + central estimate
            temp_legend(1) = {'Central Est.'};
            for l = 2:length(cflim)+1
                temp_legend(l) = {[num2str(cflim(l-1)/EVA.N_BootStrap),' quantile']};
            end
            cf_order = 1:1:(length(cflim)+1); 
            cf_order(2,1) = 0.5 * EVA.N_BootStrap; cf_order(2,2:end) = cflim; %assign the central estimate and the cf limits  
            cf_order = sortrows(cf_order.',2).'; %sort the indices for the RTP values in ascending order with the cf limits 
            temp_legend = temp_legend(cf_order(1,:)); %match the ordering of the labels with the cf order 
            temp_legend = repmat(temp_legend,1,length(EVA.title)); % create a text cell vector for the quantile values
            
            for m = 1:length(EVA.title)*(length(cflim)+1) table_row_title(m) = {[char(EVA.title(ceil(m/(length(cflim)+1)))),...
                    ' ',char(temp_legend(m))]}; end % concat the confidence labels and the directional /seasonal labels to one            
            table_row_title = circshift(table_row_title,[0,(length(EVA.title)-1)*(length(cflim)+1)]); %Make omni diectional label last rows by shifting the order of labels...
            
        else
            % only central estimate
            % for m = 1:length(EVA.title) table_row_title(m) = {[char(EVA.title(m)),' Central Est.']}; end
            table_row_title = EVA.title; % PDG: No need for central text if its the only estimate?
            table_row_title = circshift(table_row_title,[0,length(EVA.title)-1]); %Make omni diectional label last rows by shifting the order of labels...
            cf_order = 1;
            
        end
        
            m_table(filename,{EVA.name subtitle{1}},{'T_R [years]'},{[EVA.label,'  [',EVA.unit,']']},ylabs,table_row_title,...
                reshape(permute(squeeze(R_mp(:,2:end,nfit,cf_order(1,:))),[1 3 2]),[length(ylabs),length(table_row_title)]),...
                EVA.ascii,EVA.table,'Precision',EVA.ndec);        
end

XP = XP(end,:);   % return only constrained distributions

end

function[pdf] = pdf_funcs(ifit,x,g,par)

switch ifit
    
    case 2 % Truncated 2-p Weibull
        P0  = exp(-(g/par(2))^par(1));
        pdf = par(1)/par(2)*(x/par(2)).^(par(1)-1).*exp(-(x/par(2)).^par(1))/P0;
    case 3 % 2p Weibull
        pdf = par(1)/par(2)*((x-g)/par(2)).^(par(1)-1).*exp(-((x-g)/par(2)).^par(1));
    case 4 % Gumbel
        pdf = exp(-((x-par(1))/par(2))).*exp(-exp(-((x-par(1))/par(2))))/par(2);
    case 5 % Exponential
        pdf = 1/par(1).*exp(-((x-g)/par(1)));
    case 6 % Generalized Pareto
        pdf = 1/par(2)*(1+par(1)*(x-g)/par(2)).^-(1/par(1)+1);
    case 7 % Upper Bound Weibull
        P0  = exp(-((par(3)-g)/par(2)).^par(1));
        pdf = par(1)/par(2)/(1-P0)*(((par(3)-x)/par(2)).^(par(1)-1)).*(exp(-((par(3)-x)/par(2)).^par(1)));
        
end

end

function[est] = LinearRegression(x,y)
x = reshape(x,[],1);
y = reshape(y,[],1);
est = fminsearch(@fit_line, [1,0]);
    function err=fit_line(params)
        a = params(1);
        b = params(2);
        err = sum((y - (a*x+b)).^2);
    end
end

function[T_A] = T2TA(EVtype,T)
% converts from PDS T-year to Annual Max. T-year
if strcmpi(EVtype,'AMP')
    T_A = (1-exp(-1./T)).^-1;
else
    T_A = T; % already PDS series -> do nothing
end

end

function[CI] = ReturnGumbelConfInterval(T,data,mu,sigma,normz)
%
CI   = NaN(numel(T),size(data,2));
%
for j=1:size(data,2);
    nonan = ~isnan(data(:,j)); % PDG to support empty years:
    pcov  = param_cov(data(nonan,j),mu(j),sigma(j));
%     pcov  = param_cov(data(:,j),mu(j),sigma(j));
    for i=1:numel(T)
        
%         if T(i)~=1
%             zci = fzero(@(z)(exp(-exp(z + normz * sqrt(pcov(1,1)*z.^2 + 2*pcov(1,2)*z + pcov(2,2))/sigma(j)))...
%                 - 1 + 1/T(i)),log(-log(1-1/T(i))));
%         else
%             zci = fzero(@(z)(exp(-exp(z + normz * sqrt(pcov(1,1)*z.^2 + 2*pcov(1,2)*z + pcov(2,2))/sigma(j)))...
%                 - exp(-1)),exp(-1));
%         end
        % Above replaced by (2014-09-14):
        zci = fzero(@(z)(exp(-exp(z + normz * sqrt(pcov(1,1)*z.^2 + 2*pcov(1,2)*z + pcov(2,2))/sigma(j)))...
            - exp(-1/T(i))),log(1/T(i)));
        CI(i,j) = mu(j) - sigma(j)*zci;
    end
end

end

function[pcov] = param_cov(data,mu,sigma)
% inspired by evlike...
z     = (mu-data) ./ sigma;
expz  = exp(z);
L     = (z - log(sigma)) - expz;
L(z == Inf) = -Inf;
% Sum up the individual contributions, and return the negative log-likelihood.
% nlogL = -sum(L);
% Compute the negative hessian at the parameter values, and invert to get the observed information matrix.
nH11 = sum(expz);
nH12 = sum((z+1).*expz - 1);
nH22 = sum(z.*(z+2).*expz - (2.*z+1));
pcov =  (sigma.^2) * [nH11 -nH12; -nH12 nH22] / (nH11*nH22 - nH12*nH12);
end