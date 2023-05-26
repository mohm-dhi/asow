function [EVA] = m_extreme_vFLD(m_struct,options,varargin)
%
% DESCRIPTION:  Performs extreme value analysis and returns extreme value
%               structure.
%
% SYNTAX: extremedata = m_extreme(m_struct,options)
%
% REQUIRED INPUT:   m_struct    m_structure containing time series of
%                               parameter (Hm0, U10 etc.) or modedata
%
%                   options     structure containing control settings.
%                               Required fields in options are;
%
%                   .EVdist     integer defining extreme value distribution
%                               type:
%                               =2: Truncated Weibull
%                               =3: 2-p Weibull to excess
%                               =4: Gumbel
%                               =5: Exponential
%                               =6: Generalized Pareto (use with caution!)
%                               =7: Upper Bound Weibull
%
%                   .EVtype     Extreme value type;
%                               'AAP' = Average Annual Peaks
%                               'AMP' = Annual Maximum Peak
%                               'POT' = Peak-over-Threshold
%
%                   .EVcrit     vector with length = size(m_struct.data,2)
%                               if method=='AAP'; lambda (events/yr)
%                               if method=='AMP'; not required...
%                               if method=='POT'; threshold
%
%
% OPTIONAL INPUT -  optional fields in options structure - see
%                   m_DirDist_struct
%
%                   intereventtime      Inter-event time (hrs) for PDS
%                                       extraction (default 36hrs)
%
%                   intereventlevel     Inter-event level (fraction) for
%                                       PDS extraction (default 0.7)
%
%                   N_BootStrap         number of bootstraps (default 0)
%
%                   estimationmethod    Estimation method;
%                                       'ML': max likelihood (default)
%                                       'LS': least squares
%
%                   sqrfac              Perform EVA on events^(1/sqrfac)
%                                       (ie. set sqrfac=0.5 to fit EV
%                                       distribution to the square of
%                                      extremes). Default 1.0
%
%                   plotflag            Plotflag (vector)
%                                       ismember(0): Do not plot (overrides other options)
%                                       ismember(1): Plot (default)
%                                       ismember(2): include peak plots
%                                       ismember(3): include polar scatter plots
%                                       ismember(4): include lnN(Xmp) plot (m_extreme_plot_vFLD)
%                                       ismember(5): plot prod. of constrained subdists (m_extreme_plot_vFLD)
%                                       ismember(6): plot only directional
%
%
%                   EVadju              Bias correction in % (default 0)
%
%                   Hessian             Hessian confidence limits (0 or 1)
%
%                   rep_dir             representative directions (Gumbel only)
%
%                   use_PDS             use "PartialDurationSeries.cs" for
%                                       event selection:  JOAB
%
%                   optstr              Optional string for customizing file
%                                       naming when printing figures: JOAB
%
% Revision History
%
% ...
% HFH 2012-08-20:   Added output of maximum distribution (for use with m_optimize_Tr)
% HFH 2012-09-11:   Moved m_DirDist_struct to separate file
% HFH 2012-09-26:   Added dims field allowing directions to be skipped from
% HFH 2012-10-16:   Moved optional params to fields in m_struct (and defaults in m_DirDist_struct.m)
% PDG 2014-08-26:   rep_dir added as varargin in options
% PJe 2015-05-20:   Polar scatter plots of events (Plotflag 3 and directional)
% HFH 2018-10-08:   Changed plotflag to assume vector input (see above)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input reading
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=0;
if ~isempty(varargin)
    i = i+1;
    while i <= length(varargin)
        if     strcmpi(varargin{i},'invert') % PDG
            m_struct.data = -m_struct.data;
            m_struct.label = [m_struct.label '_{low}'];
            m_struct.item = [m_struct.item '_low'];
        else
            error('Unknown argument: %s',varargin{i})
        end
        i = i+1;
    end
end
%
if isempty(m_struct.data), error('Input m_structure.data was empty'); end

% Put everything into one extreme value analysis structure
EVA = m_catstruct(m_DirDist_struct,m_struct);       % create DirDist structure with default parameters and merge with input m_struct
if nargin>=2; EVA = m_catstruct(EVA,options); end   % Add options. Only the last occurence is used!
%
% Determine type of m_structure input (time series or storm modes)
if isfield(m_struct,'Modes') || isfield(m_struct,'indx')
    datatype = 1;   % modes or events
else
    datatype = 0;   % time series
end

if EVA.splitflag == 1, EVA.DIR = Z; end % PDG: add dir data for rose plot - (HFH) kan dette ikke klares med m_struct.DIR = varargin{i+1}; ?
%EVA.T = Tr; % PDG
if strcmp(EVA.item(1),'C') && length(EVA.item)==1, EVA.unit = strcat(EVA.unit,EVA.vref); end % PDG % modif FLD

%
if     EVA.EVdist == 2;   EVA.LongTermDist = 'Trunc. Weibull';
elseif EVA.EVdist == 3;   EVA.LongTermDist = '2-p Weibull';
elseif EVA.EVdist == 4;   EVA.LongTermDist = 'Gumbel';
elseif EVA.EVdist == 5;   EVA.LongTermDist = 'Exponential';
elseif EVA.EVdist == 6;   EVA.LongTermDist = 'Generalized Pareto';
elseif EVA.EVdist == 7;   EVA.LongTermDist = 'Weibull Upper Bound';
else                    error('Unknown Long-Term Distribution');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extreme event extraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if datatype==0
    % extract events from time series input data type
    events = m_PDS(m_struct,EVA.intereventtime,EVA.intereventlevel,EVA.use_PDS);
    % polar plot of all events
    if ismember(3,EVA.plotflag) && isfield(EVA,'directional')
        if EVA.directional
            temp = events;
            temp.title{1} = 'All Omni events';
            temp.data = temp.data(:,1);
            Y.time = temp.time;
            Y.label = temp.sublabel;
            Y.data = temp.subdata(temp.indx(:,1));
            Y.item = temp.sublabel;
            Y.legend = temp.legend;
            Y.isdir = 1;
            Y.group = {'all'};
            Y.bins = temp.subbins;
            
            m_polarscatter(temp,Y,'density','dirbins')
            
            temp.title{1} = 'All directional events';
            temp.data = events.data(:,2:end);
            nc = size(temp.data,2);
            Y.data =[];
            for i = 1:nc
                nr(i) = size(temp.data,1);
                for j = 1:nr(i)
                    if isnan(temp.indx(j,i+1))
                        Y.data(j,i) = nan;
                    else
                        Y.data(j,i) =temp.subdata(temp.indx(j,i+1));
                    end
                end
            end
            Y.time = temp.time;
            Y.data = Y.data(:);
            temp.data = temp.data(:);
            m_polarscatter(temp,Y,'density','dirbins')
            
        end
    end
    
else
    % use modes or events from m_struct directly
    events.time = m_struct.time;
    events.data = m_struct.data;
    events.Nyrs = m_struct.Nyrs;
    if isfield(m_struct,'indx')
        events.indx = m_struct.indx;
    end
    
    % add some additional fields specific for mode data
    if isfield(m_struct,'Modes')
        % put modes and Neqv in EVA structure for subsequent convolution
        tmp = [m_struct.Modes(:,1).Mode];
        EVA.mode_N_shape = [tmp(~isnan(tmp));[m_struct.Modes(~isnan(tmp),1).N_eqv]]';
        % add shape parameter if present (this is the case for Forristall_C,
        % where it has been computed for each storm sea state and is stored at
        % the storm peak)
        if strcmp(m_struct.ShortTermDist,'Forristall_C')
            EVA.mode_N_shape = [EVA.mode_N_shape,[m_struct.Modes(~isnan(tmp),1).beta_peak]'];
        end
        if strcmp(m_struct.ShortTermDist,'Glukhovskiy_H')
            EVA.mode_N_shape = [EVA.mode_N_shape,[m_struct.Modes(~isnan(tmp),1).kappa_peak]'];
        end
        % also put short-term distribution identifier in EVA structure
        EVA.ShortTermDist = m_struct.ShortTermDist;
    end
end
if ~isfield(events,'indx')
    events.indx = nan(size(events.data));
end
events.data(isnan(events.data)) = 0;        % replace NaNs with zeros (for sorting)
events.data = circshift(events.data,[0,-1]);% shift omni-values to last column (this is how DirDist wants data)
events.indx = circshift(events.indx,[0,-1]);
EVA.data = circshift(EVA.data,[0,-1]);% PDG
nc = size(events.data,2);                   % number of columns in data (ndir + 1)

% Perform directional fitting
EVA.Nyrs      = events.Nyrs;
% EVA.EVadju    = EVadju; % PDG
% EVA.Tcon      = Tcon; % PDG
% EVA.intereventtime = intereventtime; % PDG
% EVA.intereventlevel = intereventlevel; % PDG
% EVA.splitflag = splitflag; % PDG

% % Set some optional inputs
% if exist('EstimationMethod','var'); EVA.method = EstimationMethod;   end
% if exist('N_BootStrap','var');      EVA.N_BootStrap = N_BootStrap;   end
% if exist('sqrfac','var');           EVA.sqrfac = sqrfac;             end
% if exist('skipdirfrac','var');      EVA.skipdirfrac = skipdirfrac;   end
% if exist('constfac','var');         EVA.constfac = constfac;         end
% if exist('plotflag','var');         EVA.plotflag = plotflag;         end
if ~isfield(EVA,'dims');            EVA.dims = 1:size(EVA.data,2);   end    % allow user to specify directions to include in analysis

%if exist('ConfLimits','var');       EVA.ConfLimits = ConfLimits;     end

% reselection from available events based on selected method
% we use the same lambda for all directions (perhaps this should also be randomly sampled from Poisson process)
%if length(EVA.EVcrit)==1; EVA.EVcrit = repmat(EVA.EVcrit,1,nc); end % expand lambda to number of columns
if length(EVA.EVcrit)>1; disp('WARNING - m_extreme: Different lambdas for the different directions does not work at the moment...'); end
EVA.EVcrit = repmat(EVA.EVcrit(1),1,nc);
N = round(events.Nyrs*EVA.EVcrit);                          % number of events to extract for each column
for c = 1:nc
    [X_all(:,c),iAll] = sort(events.data(:,c),'descend');                    % sort all events descending
    X_all_time(:,c) = events.time(iAll);
    X_all_i(:,c) = events.indx(iAll,c);
end
switch EVA.EVtype
    
    case 'AAP'
        % This method uses a fixed average number of exceedances per year to
        % define threshold. The same number of events are extracted for the
        % bootstrap resampling. Alternatively one could consider sampling the
        % number of exceedances from the corresponding Poisson process
        
        X = NaN(max(round(events.Nyrs*EVA.EVcrit)),nc);         % initialize array of selected events
        Xtime = X.*nan;
        Xindx = X.*nan;
        for c=1:nc
            X(1:N(c),c) = X_all(1:N(c),c);
            if datatype==0
                Xtime(1:N(c),c) =  X_all_time(1:N(c),c);
                Xindx(1:N(c),c) =  X_all_i(1:N(c),c);
            end
        end          % load selected events into array
        
        i_end   = (size(X,1)-sum(isnan(X)));                    % find indices of last numeric element for each direction
        istp = 4; g = X(sub2ind(size(X),i_end,1:size(X,2)))-... % use 'slope' of events calculated for lowest istp points to set appropriate threshold (ML fits quite sensitive to this)
            (X(sub2ind(size(X),i_end-istp,1:size(X,2)))-X(sub2ind(size(X),i_end,1:size(X,2))))./istp/2;
        
    case 'POT' % PDG
        
        if EVA.splitflag ~= 0; error('EVtype ''POT'' implemented for omni/all conditions only!'); end % Since each column gets different length or NaNs value
        
        X = NaN(sum(X_all(:,end)>EVA.EVcrit(end)),nc);          % initialize array of selected events (based on omni)
        Xtime = X.*nan;
        Xindx = X.*nan;
        for c=1:nc
            X(1:sum(X_all(:,c)>EVA.EVcrit(c)),c) = X_all(find(X_all(:,c)>EVA.EVcrit(c)),c);
            if datatype==0
                Xtime(1:sum(X_all(:,c)>EVA.EVcrit(c)),c) = X_all_time(find(X_all(:,c)>EVA.EVcrit(c)),c);
                Xindx(1:sum(X_all(:,c)>EVA.EVcrit(c)),c) = X_all_i(find(X_all(:,c)>EVA.EVcrit(c)),c);
            end
        end  % load selected events into array
        
        g = ones(1,nc) * EVA.EVcrit(1);  %!!!!!!!!!!!!!!!! skal ændres hvis dette udvides til at håndtere multidim.
        
    case 'AMP'
        if EVA.EVdist ~= 4; warning('EVdist should be Gumbel(4) for use with EVtype ''AMP''!'); end % Pt...? PDG
        X = NaN(round(EVA.Nyrs),nc);     % initialize array of selected events (based on omni)
        Xtime = X.*nan;
        Xindx = X.*nan;
        V = datevec(events.time);
        n = size(V,1);
        EVA.startdate = [V(1,1),EVA.startdate];
        %
        % make pseudo dates for events by shifting data prior to start date
        % to the end and then move everything backwards in time to the Jan.
        % 1 of the first ýear.
        i0 = find(events.time<datenum(EVA.startdate),1,'last'); if isempty(i0); i0=0; end
        pseudo_time = datevec(datenum([[V(i0+1:n,1);V(1:i0,1)+EVA.Nyrs],[V(i0+1:n,2:6);V(1:i0,2:6)]]) ...
            - (datenum(EVA.startdate)-datenum([EVA.startdate(1),1,1])));
        pseudo_data = events.data([i0+1:n,1:i0],:);
        pseudo_indx = events.indx([i0+1:n,1:i0],:);
        %
        % loop trough years and extract maximum in each column for
        % particular year
        
        for i=1:round(EVA.Nyrs)
            %             X(i,:) = max(pseudo_data(pseudo_time(:,1) == i-1+EVA.startdate(1),:)); % PDG: breaks down for empty years!
            I0 = pseudo_data(pseudo_time(:,1) == i-1+EVA.startdate(1),:);
            indx0 = pseudo_indx(pseudo_time(:,1) == i-1+EVA.startdate(1),:);
            T0 = pseudo_time(pseudo_time(:,1) == i-1+EVA.startdate(1),:); % PDG
            if isempty(I0)
                X(i,:) = NaN;
                Xindx(i,:) = NaN;
            else
                for jj= 1:nc
                    if contains(EVA.item,'low')
                        [X(i,jj),Ipt] = max(I0(find(I0(:,jj) ~= 0),jj)); % find maximum non-zero "negative" number
                        Ipt = find(I0(:,jj) == X(i,jj)); % nan's messing with desired command
                    else
                        [X(i,jj),Ipt] = max(I0(:,jj));
                    end
                    if X(i,jj) ==0
                        continue
                    end
                    %                     if isfinite(sum(Xindx)) % PDG: quick add to curcumvie error if Xindx are nans (not i structure)
                    if datatype==0
                        Xindx(i,jj) = indx0(Ipt,jj);
                        Xtime(i,jj) = EVA.time(Xindx(i,jj));
                    end
                end
            end
            if ~isempty(T0)
                EVA.Tx(i,:) = datenum(T0(Ipt,:)) + (datenum(EVA.startdate)-datenum([EVA.startdate(1),1,1])); % PDG
            else
                EVA.Tx(i,:) = NaN;
            end
        end
        
        % Find represented directional sectors by accepting that only up to 10% of years have no data in a sector, and where 10% of the data
        % points are more than 5% of the max value and make reduced temporary structure containing only relevant directions
        if EVA.rep_dir
            if size(X,2) == 1 % for WL
                EVA.rep_dirs = 1;
            else
                EVA.rep_dirs = find(sum(i0==0)<0.1*EVA.Nyrs & sum(X<=0.05*max(max(X)))<0.1*EVA.Nyrs); % rep_dir - PDG
                X       = X(:,EVA.rep_dirs);
                EVA.dims = 1:length(EVA.rep_dirs);
                EVA.EVcrit = EVA.EVcrit(EVA.rep_dirs);
                nc = length(EVA.rep_dirs);
            end
        end
        X(isnan(X)) = 0; % PDG: Fix if there is no data for 1 particular year - but really one should not use AMP then!! %X(X==0) = min(min(X(X>0)));
        if sum(sum(X==0)) >0, warning('Likely missing storm(s) (rerun modes with lower threshold)'), end % PDG warning for empty modes
        %
        g = zeros(1,nc);               % Set location g = 0 for AMP method
        
        clear pseudo_time pseudo_data i0 V n
        
end
if ~isfield(m_struct,'Modes')
    EVA.Xtime = Xtime;
    EVA.Xindx = Xindx;
end

% polar scatter plot of selected event data added by PJe
if ismember(3,EVA.plotflag) && isfield(EVA,'directional') && ~isfield(m_struct,'Modes')
    if EVA.directional
        temp = events;
        temp.title{1} = 'Directional extreme events';
        temp.data = X(:,1:end-1);
        nc = size(X,2)-1;
        clear Y;
        Y.data =nan(size(EVA.Xindx,1),nc);
        for i = 1:nc
            Y.data(~isnan(EVA.Xindx(:,i)),i) =  EVA.dirstruc.data(EVA.Xindx(~isnan(EVA.Xindx(:,i)),i),1);
        end
        Y.time = temp.time;
        Y.data = Y.data(:);
        temp.data = temp.data(:);
        Y.label = temp.sublabel;
        Y.item = temp.sublabel;
        Y.legend = temp.legend;
        Y.isdir = 1;
        Y.group = {'all'};
        Y.bins = temp.subbins;
        
        %m_polarscatter(temp,Y)
        m_polarscatter(temp,Y,'dirbins')
        
        temp.data = X(:,end);
        
        temp.title{1} = 'Omni extreme events';
        Y.data =  EVA.dirstruc.data(EVA.Xindx(:,end),1);
        Y.data = Y.data(:);
        m_polarscatter(temp,Y,'dirbins')
    end
end

clear events
[~,EVA.Location]  = ndgrid(1,g,1:EVA.N_BootStrap+1);            % and expand fixed treshold to n_bootstrp+1 dimensions
[~,EVA.lambda]    = ndgrid(1,EVA.EVcrit,1:EVA.N_BootStrap+1);   % also expand to n_bootstrp+1 dimensions

% expand array of events with randomly sampled (with replacement) from original data
X_bootstrap = X;                                        % orig. values in first layer
X_bootstrap(:,:,2:EVA.N_BootStrap+1) = NaN;             % init. array of bootstrapped values
cSub = meshgrid(1:size(X,2),1:size(X,1)); cSub=cSub(:); % column sub-indices
for i=1:EVA.N_BootStrap
    rSub  = randi(size(X,1),size(cSub));                % random row indices
    X_bootstrap(:,:,i+1) = reshape(X(sub2ind(size(X),rSub,cSub)),size(X)); % add events from randomly sampled indices
end
EVA.X      = X_bootstrap;                         % use bootstrap events as input

% make plot of events used in EVA
if ismember(2,EVA.plotflag) && ~datatype; m_peak_plot(EVA); end
%
% Apply peak adjustment PDG
EVA.X        = EVA.X*(1+EVA.EVadju/100);
EVA.Location = EVA.Location*(1+EVA.EVadju/100);
%
% raise to the power of 1/sqrfac
EVA.X        = EVA.X.^(1/EVA.sqrfac);
EVA.Location = EVA.Location.^(1/EVA.sqrfac);
%
% call DirDist function to fit distributions and constrain to omni if required
%--------------------------------------------------------------------------
% JOR: Rare cases where not enough data fill the vals/per yer matrix causes
% zeros to ruin the fit...
[rr,cc] = find(EVA.X==0);
   if all(size(EVA.data,2) == 13 && size(cc,2)<=5) % JOGR monthly minimum: 5 exceptions for 1 event year
    disp('WARNING - Dirs with less events per year NOT being neglected...')
    elseif ~isempty(cc)  % PDG added if loop to avoid disp below when not relevant
%    if all(size(EVA.data,2) == 13 && size(cc,2)>=5) % JOGR monthly minimum exception
    cc2 = unique(cc); % PDG: not sure I follow this...
    EVA.X(:,cc2) = 0;
    EVA.Location(1,cc2) = 0;
    disp('WARNING - Dirs with less events per year being neglected...')

end, clear cc
%--------------------------------------------------------------------------

EVA          = m_DirDist(EVA);

% Make probability distribution plot and calculate return period values
% [EVA]        = m_extreme_plot_vFLD(EVA); % PDG
[EVA]        = m_extreme_plot_HL09(EVA); % PDG
%
% EVA.R_mp(:,:,2,1)'
% EVA.R_mx(:,:,2,1)'

% rep_dir (reorder into full matrix for output - relevant fields only)
EVA.R_mp(:,1,:) = []; % Remove first col
EVA.R_mx(:,1,:) = []; % Remove first col
if EVA.rep_dir && length(EVA.rep_dirs)<size(EVA.data,2)
    not_rep_dir = setdiff(1:size(EVA.data,2),EVA.rep_dirs);
    
    EVA.X(:,EVA.rep_dirs) = EVA.X;
    EVA.X(:,not_rep_dir)  = nan(size(EVA.X,1),length(not_rep_dir));
    
    EVA.UnconstrainedParams(:,EVA.rep_dirs) = EVA.UnconstrainedParams;
    EVA.UnconstrainedParams(:,not_rep_dir)  = nan(size(EVA.UnconstrainedParams,1),length(not_rep_dir));
    
    if ~isempty(EVA.ConstrainedParams)
        EVA.ConstrainedParams(:,EVA.rep_dirs) = EVA.ConstrainedParams;
        EVA.ConstrainedParams(:,not_rep_dir)  = nan(size(EVA.ConstrainedParams,1),length(not_rep_dir));
    end
    
    EVA.R_mp(:,EVA.rep_dirs,:) = EVA.R_mp(:,:,:);
    EVA.R_mp(:,not_rep_dir,:) = nan(length(EVA.T),length(not_rep_dir),2);
    
    EVA.R_mx(:,EVA.rep_dirs,:) = EVA.R_mx(:,:,:);
    EVA.R_mx(:,not_rep_dir,:) = nan(length(EVA.T),length(not_rep_dir),2);
    
    if EVA.Hessian && EVA.EVdist == 4 % (Gumbel only)
        EVA.R_LB(:,EVA.rep_dirs) = EVA.R_LB;
        EVA.R_LB(:,not_rep_dir)  = nan(size(EVA.R_LB,1),length(not_rep_dir));
        
        EVA.R_UB(:,EVA.rep_dirs) = EVA.R_UB;
        EVA.R_UB(:,not_rep_dir)  = nan(size(EVA.R_UB,1),length(not_rep_dir));
    end
    
end

end


% HFH: do we need this? PDG: Yes
function m_peak_plot(EVA) % PDG

EVA.title = circshift(EVA.title,[0,-1]); % only reorder once

nc = size(EVA.X,2);
for i = 1:nc % PDG: only for omni?
    
    nr(i) = size(EVA.X,1);
    %     for j = 1:nr(i)
    %         It(j) = find(EVA.data(:,i)==EVA.X(j,i)); % PDG: Very dirty - and fails if find more than 1 value (low precision) ()but opnlu for plotting though
    %     end
    It = EVA.Xindx(:,i);
    Ifin = find(isfinite(It));
    
    
    % Filename
    filename = strcat(EVA.name,'_Extreme_',EVA.item,'_',EVA.legend,EVA.ttt_str,EVA.optstr);
    filename(isspace(filename)) = '_';
    
    % Define string
    if     strcmp(EVA.EVtype,'POT'), str = 'Threshold value';
    elseif strcmp(EVA.EVtype,'AMP'), str = ['Annual max peaks'];
    elseif strcmp(EVA.EVtype,'AAP'), str = 'Avg. annual peaks';
    end
    
    %flip inverted data for plotting
    if contains(EVA.item,'low') && contains(EVA.item,'Tair')
        EVA.data = -EVA.data;
        EVA.X = -EVA.X;
        title_str = 'Peaks_{low}';
        bin_diff = EVA.bins(2)-EVA.bins(1);
        bin_1 = EVA.bins(1); EVA.bins = bin_1-bin_diff*2:bin_diff:EVA.bins(end);
    else
        title_str = 'Peaks';
    end
    
    % Plot peak timeseries
    plot(EVA.time,EVA.data(:,i),'color',EVA.ColorOrder(2,:),'Marker','.','LineStyle','none'), hold on % ,'LineWidth',1.5 ,'MarkerSize',3
    plot(EVA.time(It(Ifin)),EVA.X(Ifin,i),'color',EVA.ColorOrder(1,:),'LineStyle','none','Marker','o','MarkerSize',3,'LineWidth',1.5), grid on
    %     plot(EVA.Tx,EVA.X(:,i),'color',EVA.ColorOrder(1,:),'LineStyle','none','Marker','o','MarkerSize',3,'LineWidth',1.5), grid on
    title([strcat(EVA.name,EVA.xyz_str),strcat(title_str,EVA.ttt_str,{' '},EVA.legend,{' - '},EVA.title{i})],'FontWeight','normal')
    %     legend(['Data' ' (N = ' num2str(length(EVA.time)) ', N_{years} = ' num2str(EVA.Nyrs,'%2.1f') ')'],['Peaks' ' (N = ' num2str(nr(i)) ', \lambda = ' num2str(EVA.lambda(i),'%02.1f') ', g = ' num2str(EVA.Location(i),'%02.2f') ')'],'location','SouthEast')
    legend(['Data' ' (N = ' num2str(length(EVA.time)) ', N_{years} = ' num2str(EVA.Nyrs,'%2.1f') ')'],['Peaks' ' (N = ' num2str(nr(i)) ', \lambda = ' num2str(EVA.lambda(i),'%02.1f') ', g = ' num2str(EVA.Location(end),'%02.2f') ')'],'location','SouthEast')
    set(gca,'ylim',[EVA.bins(1) EVA.bins(end)]), set(gca,'ytick',EVA.bins), ylabel(strcat(EVA.label,' [',EVA.unit,']'))
    [xticks,xlabels] = m_date_tick([EVA.ttt(1) EVA.ttt(2)]); set(gca,'xlim',[EVA.ttt(1) EVA.ttt(2)]), m_xticklabel_rotate(xticks,45,xlabels)
    if strcmpi(EVA.item,'WL_low')
        set(gca,'YTickLabel',flipud(get(gca,'YTickLabel'))),
    end % flip YTickLabel for low water level
    
    % Size and Position
    fpos = get(gcf,'Position'); set(gcf,'Position',[fpos(1)/2   fpos(2)/2   2*fpos(3)   fpos(4)],'PaperPositionMode','auto')
    set(gca,'Position',[0.07 0.17 0.9 0.73])
    ylabel_pos = get(get(gca,'ylabel'),'Position'); set(get(gca,'ylabel'),'position',[0.45*ylabel_pos(1) ylabel_pos(2) ylabel_pos(3)]);
    xlabel_pos = get(get(gca,'xlabel'),'Position'); set(get(gca,'xlabel'),'position',[xlabel_pos(1) 0.7*xlabel_pos(2) xlabel_pos(3)]);
    
    % Annotate and print
    x_pos1 = EVA.ttt(1)+0.012*(EVA.ttt(2)-EVA.ttt(1));
    if strcmp(EVA.EVtype,'AMP')
        rectangle('Position',[EVA.ttt(1)+0.004*(EVA.ttt(2)-EVA.ttt(1)) EVA.bins(1)+0.012*(EVA.bins(end)-EVA.bins(1)) 0.18*(EVA.ttt(2)-EVA.ttt(1)) 0.29*(EVA.bins(end)-EVA.bins(1))],'FaceColor','w')
        text(x_pos1,EVA.bins(1)+0.26*(EVA.bins(end)-EVA.bins(1)),['Start date (mm-dd) = ' datestr(datenum([2000 EVA.startdate(2:3)]),'mm-dd')])
    else
        rectangle('Position',[EVA.ttt(1)+0.004*(EVA.ttt(2)-EVA.ttt(1)) EVA.bins(1)+0.012*(EVA.bins(end)-EVA.bins(1)) 0.18*(EVA.ttt(2)-EVA.ttt(1)) 0.22*(EVA.bins(end)-EVA.bins(1))],'FaceColor','w')
    end
    if i == nc, num = 0; else, num = i; end % omni is last!
    text(x_pos1,EVA.bins(1)+0.19*(EVA.bins(end)-EVA.bins(1)),[str ' = ' num2str(EVA.EVcrit(end),'%02.2f')])
    text(x_pos1,EVA.bins(1)+0.12*(EVA.bins(end)-EVA.bins(1)),['Inter-event time = ' num2str(EVA.intereventtime,'%2.1f') 'h'])
    text(x_pos1,EVA.bins(1)+0.05*(EVA.bins(end)-EVA.bins(1)),['Inter-event level = ' num2str(EVA.intereventlevel,'%2.1f')])
    substr = [EVA.EVtype '_' num2str(EVA.EVcrit(end)) '_' num2str(EVA.intereventtime,'%2.1f') 'h' '_' num2str(EVA.intereventlevel,'%2.1f')];
    fname = strcat(filename,'_',substr,'_',num2str(num,'%02.0f'),'_',EVA.title{i},'_Time',EVA.figure);
    fname = m_filename(fname);
    print(gcf,'-dpng',EVA.reso,fname)
    close()
    
    % return to previous values for rest of script to work as expected
    if contains(EVA.item,'low') && contains(EVA.item,'Tair')
        EVA.data = -EVA.data;
        EVA.X = -EVA.X;
        title_str = 'Peaks';
        EVA.bins = bin_1:bin_diff:EVA.bins(end);
    else
        title_str = 'Peaks';
    end
    
    % Plot peak rose
    if EVA.subseries && strcmp(EVA.subtype,'directional') && i==nc
        % m_rose_plot(EVA.subdata(It(Ifin)),EVA.X(Ifin,i),'n',length(EVA.subbins)-1,'di',EVA.bins(EVA.bins>EVA.bins(end)/4),'lablegend',[{strcat(EVA.label,' (',EVA.unit,')')},{strcat(' vs. ',EVA.sublabel)}],'labtitle',[strcat(EVA.name,EVA.xyz_str),strcat('Peak Rose Plot',EVA.ttt_str,{' '},EVA.legend,{' - '},EVA.title{i})],'dtype','non-meteo');
        if EVA.subbins(1)~=0 && size(EVA.subbins,2)<=3 % JOAB edited for two-directional bins non-equidistantly spaced.
            m_rose_plot(EVA.subdata(It(Ifin)),EVA.X(Ifin,i),'n',length(EVA.subbins)-1,'di',EVA.bins(EVA.bins>EVA.bins(end)/4),'lablegend',[{strcat(EVA.label,' [',EVA.unit,']')},{strcat(' vs. ',EVA.sublabel)}],'labtitle',[strcat(EVA.name,EVA.xyz_str),strcat('Peak Rose Plot',EVA.ttt_str,{' '},EVA.legend,{' - '},EVA.title{i})],'dtype','standard','ci',[],'cmap',EVA.ColorMap,'Ay',[EVA.subbins - 0.5*(EVA.subbins(2)-EVA.subbins(1))-0.5*EVA.subbins(1)]);
        else
            m_rose_plot(EVA.subdata(It(Ifin)),EVA.X(Ifin,i),'n',length(EVA.subbins)-1,'di',EVA.bins(EVA.bins>EVA.bins(end)/4),'lablegend',[{strcat(EVA.label,' [',EVA.unit,']')},{strcat(' vs. ',EVA.sublabel)}],'labtitle',[strcat(EVA.name,EVA.xyz_str),strcat('Peak Rose Plot',EVA.ttt_str,{' '},EVA.legend,{' - '},EVA.title{i})],'dtype','standard','ci',[],'cmap',EVA.ColorMap,'Ay',[EVA.subbins - 0.5*(EVA.subbins(2)-EVA.subbins(1))]);
        end
        % m_rose_plot(X.data(nonans,i),Y.data(nonans,i),'n',c,'di',[Y.bins(find(Y.bins>=Ymin,1,'first')) Y.bins(Y.bins>=Ymin & Y.bins<=Y.bins(end)) Y.bins(end)],'lablegend',[X.legend,strcat({'N = '},int2str(N)),strcat(X.label,' [',X.unit,']'),strcat(Y.label,' [',Y.unit,']')],'labtitle',[heading_rose(1) strcat(heading_rose(2),{' '},X.title{i})],'dtype','standard','ci',ci,'cmap',X.ColorMap,'Ay',[X.bins - 0.5*(X.bins(2)-X.bins(1))]);
        print(gcf,'-dpng',EVA.reso,strcat(filename,'_',substr,'_',num2str(num,'%02.0f'),'_',EVA.title{i},'_Rose',EVA.figure)), close
        
        % PDG: added directional 'rose-plot' - showing highest events for each DIR sector:
        m_rose_plot(EVA.subdata(reshape(EVA.Xindx(:,1:end-1),[size(EVA.Xindx,1)*(size(EVA.Xindx,2)-1),1])),reshape(EVA.X(:,1:end-1),[size(EVA.X,1)*(size(EVA.X,2)-1),1]),'n',length(EVA.subbins)-1,'di',EVA.bins(EVA.bins>EVA.bins(end)/4),'lablegend',[{strcat(EVA.label,' [',EVA.unit,']')},{strcat(' vs. ',EVA.sublabel)}],'labtitle',[strcat(EVA.name,EVA.xyz_str),strcat('Peak Rose Plot',EVA.ttt_str,{' '},EVA.legend,{' - '},'Directional')],'dtype','standard','ci',[100/(length(EVA.subbins)-1)],'cmap',EVA.ColorMap,'Ay',[EVA.subbins - 0.5*(EVA.subbins(2)-EVA.subbins(1))]);
        print(gcf,'-dpng',EVA.reso,strcat(filename,'_',substr,'_',num2str(num,'%02.0f'),'_',EVA.title{i},'_Rose_Directional',EVA.figure)), close
        
    end
end

end
