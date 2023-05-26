function [modedata] = m_mode_vFLD(varargin)

% DESCRIPTION: Find storm wave modes based on peak events and return
%              modedata structure.
%
% SYNTAX:           modedata = m_mode(Hm0,T0,varargin)
%
% REQUIRED INPUT:   Hm0,T02  (m_structures)
%
% OPTIONAL INPUT:   m_structures of other environmental variables followed
%                   by pairwise 'string',variable inputs. The possible
%                   optional inputs can be found in the list of default
%                   parameters in m_mode.m
%
% EXAMPLE:  modedata = m_mode(Hm0,T02,Tp,MWD,WL,'intereventtime',18,...
%           'output_types',1:2,'shorttermdist','Forristall_H')
%
% Copyright DHI ©, Denmark, All rights reserved
% 2012-06-15 - Hans Fabricius Hansen, hfh@dhigroup.com
%
% Revision History
% HFH 2012-06-15:	First release
% PDG 2012-06-30:   Search for PDG below
% HFH 2012-10-16:   Added 'rounded2nearest' field. Also removed use of Hm0.ttt(3)
% HFH 2014-05-06:   Misc. fixes for depth-dependent distributions
% FLD/
% HFH 2016-03-14:   Addition of 'Gauss Bell' equivalent storm model
%                   including associated (weighted mean) parameters.
%                   Introduction of matlab-based m_decluster which limits
%                   storm time series to up-crossing of StormCutOff level.
% HFH 2016-04-28:   Added optional interpolation across data gaps
% HFH 2016-09-07:   More generic input of associated variables
% HFH 2018-11-13:   Major reworking to make it compatible with J-EVA
% HFH 2018-11-23:   Increased MaxIter for fminsearch
% HFH 2018-12-05:   Additional check of time syncronization
% HFH 2018-12-20:   Fixed wrong ln(sigma) start_guess for FitSTDStorm
% JOAB 2019-25-11:  Fixed a bug with shorttermdist naming not working
% properly for Forristall_H wave height distributions
% Default parameters
shorttermdist   = 'Rayleigh_H';
intereventtime  = 36;   % inter-event time in hours
intereventlevel = 0.75; % inter-event level (required drop between adjacent peaks in fraction of lowest peak value)
StormCutOff     = [0.75,0.75];  % include time steps from first up-crossing of StormCutOff(1) to last downcrossing of StormCutOff(2) into TSvalues
plotflag        = 0;    % plotflag=1 triggers figure plotting and saving (slow)
dur             = [];   % force the sea state duration to take this value (in seconds)
OutputFolder    = '.\';   % specify an output folder (default will be matlab working directory)
LowerLimitPrct  = [];   % can be used to remove events less than LowerLimitPrct times max. Hm0 peak.
fontsize        = 12;   % fontsize used on figures
rho             = [];   % parameter required for Næss wave height distribution
slope           = [];   % seabed slope required for Battjes and Groenendijk wave heights
EVadju          = 0; % PDG: peak adjustment in prct. as Hm0.data=(Hm0.data*(1+EVadju/100);
return_only_maximum_mode = 0; % return_only_maximum_mode=n returns n largest events (by peak Hm0). Ignored if return_only_maximum_mode=0
min_Hm0_peak    = 0;   	% only return events with peak Hm0 above this value
interpol_hrs    = 0;    % close gaps less than interpol_hrs by interpolation
interpol_mthd   = 'linear'; % use linear interpolation to close gaps (see interp1.m for other options)
remove_TS       = 0;    % remove TSvalues field in order to reduce file size
StartGuess      = [];
mnsz_lnlnP      = 1;    % 0: minimize p. 1:  minimize difference in ln(-ln(P))
dx              = 0.2;  % spacing on X vector for computing CDF
output_types    = 1;    % ismember(1,output_types): output for m_extreme
                        % ismember(2,output_types): output for J-EVA
Forristall_C_dim = 3; % use Forristall 3D crest height distribution by default
use_PDS         = 0; % option to use old PartialDurationSeries tool for declustering (not recommended)
decluster_item  = 'Hm0'; % option to specify other item than Hm0 for declustering (eg. wind speed, current speed)
time_tol        = 60;   % maximum tolerance on time syncronization (in seconds)
%
%% find m_structure inputs
cc=1;
in_mstr = false(1,length(varargin));
for i=1:length(varargin)
    if isstruct(varargin{i}) && isfield(varargin{i},'item')
        in_mstr(i) = true;
        k_item = find([isspace(varargin{i}.item) length(varargin{i}.item)],1,'first')-1; % PDG: to support items with 'sub-name' (space etc.)
        varargin{i}.item = varargin{i}.item(1:k_item);
        fprintf('reading %s\n',varargin{i}.item);
        eval([varargin{i}.item,'=varargin{i};'])
        itemNames{cc}  = varargin{i}.item;
        labels{cc}     = varargin{i}.label;
        units{cc}      = varargin{i}.unit;
        if strfind(varargin{i}.unit,'\circN')
            itemDir(cc)=1;
        else
            itemDir(cc)=0;
        end
        cc=cc+1;
    elseif strcmp(varargin{i},'IEL')
        intereventlevel = varargin{i+1}; cc=cc+1;
    elseif strcmp(varargin{i},'IET')
        intereventtime = varargin{i+1}; cc=cc+1; 											
    end
end
%
assign_argin(varargin(~in_mstr))  % assign optional 'string',var inputs
parm = upper(shorttermdist(end)); % either 'H' or 'C'
if ~isempty(OutputFolder), if exist(OutputFolder,'dir') ~= 7, mkdir(OutputFolder); OutputFolder = strcat(OutputFolder,'\'); end; end % PDG
Hm0.data = Hm0.data*(1+EVadju/100); % Apply peak adjustment PDG
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Input reading and checking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if length(unique(Hm0.group(2:end)))>1
    error('Do not mix directional and monthly subsets');
end
%
% Mean wave period definition
if size(T02.data,1)~=size(Hm0.data,1); error('Size of T02 does not match size of Hm0'); end
if isempty(dur); dur = median(Hm0.time(2:end)-Hm0.time(1:end-1))*24*3600;end % sea state duration in seconds % PDG
disp(['m_mode: Sea state duration is ',int2str(dur/60),' minutes']);
disp(['m_mode: Water depth is ',num2str(abs(Hm0.xyz(3)),'%5.1f'),' meters']); % PDG 'wd' replaced with abs(Hm0.xyz(3))
%
% add additional parameters if defined
for i=1:length(itemNames)
    eval(['m_str = ',itemNames{i},';']);
    if size(m_str.data,1)~=size(Hm0.data,1); error('Size of %s does not match size of Hm0',itemNames{i}); end
    if any(abs(m_str.time-Hm0.time)>time_tol/24/3600); error('Time steps of %s deviates by more than %i secs from time steps of Hm0',itemNames{i},time_tol); end
    itemInfo.(itemNames{i}) = rmfield(m_str,{'data','time','MarkerOrder','LineOrder'}); %'dfs_name','dfs_unit','dfs_type'
end

V = datevec(Hm0.time); yrs = unique(V(:,1));
if min(Hm0.time) >= datenum([min(yrs) 03 01]), yrs = yrs(2:end  ); end
if max(Hm0.time) <= datenum([max(yrs) 02 29]), yrs = yrs(1:end-1); end
Nyrs = (Hm0.time(end)-Hm0.time(1)-length(yrs(find(eomday(yrs,2)==29))))/365;
%if isempty(Npeaks); Npeaks = round(4*Nyrs); end % default average annual number of events = 4

% PDG: Filename and subtitle
subtitle = strcat({' ('},datestr(Hm0.ttt(1),29),{' - '},datestr(Hm0.ttt(2),29),{') '});
filename = strcat(Hm0.name,'_Mode_',Hm0.legend,subtitle{1},'_IET=',num2str(intereventtime,'%2.1f'),'h_','IEL=',num2str(intereventlevel,'%2.2f'),'_',shorttermdist);
filename(isspace(filename)) = '_';
if strcmp(parm,'C') || strcmp(parm,'T')
    if exist('WL','var')
        filename = strcat(filename,'_',WL.vref);
    else
        filename = strcat(filename,'_','SWL');
    end
end
%
% Forristall crest height distribution
if strcmp(shorttermdist,'Forristall_C')
    DistIn.dim = Forristall_C_dim;
    DistIn.alfa = [];
    DistIn.beta = [];
end
%
% Forristall Wave-height distribution
if strcmp(shorttermdist,'Forristall_H')
%     DistIn.dim = Forristall_H;
    DistIn.alfa = [];
    DistIn.beta = [];
end
%
% Næss Wave Height Distribution
if strcmp(shorttermdist,'Naess_H')
    if isempty(rho); error('parameter rho must be defined when Næss wave height distribution is used'); end
    DistIn.rho = rho;
    filename = strcat(filename,'_rho=',num2str(rho));
end
%
% Battjes & Groenendijk Wave Height Distribution
if strcmp(shorttermdist,'Battjes_Groenendijk_H')
    if isempty(slope); error('parameter slope must be defined when Battjes & Groenendijk wave height distribution is used'); end
    DistIn.slope = slope;
    filename = strcat(filename,'_slope=',num2str(slope));
    DistIn.d     = abs(Hm0.xyz(3));
end
if size(Hm0.data,2)>1
    filename = strcat(filename,['_',char(Hm0.group{2})]);
end
%
% Weibull Wave Height Distribution
if strcmp(shorttermdist,'Weibull_H')
    DistIn.alpha = alpha;
    DistIn.beta = beta;
    filename = strcat(filename,'_alpha=',num2str(alpha),'_beta=',num2str(beta));
end
%
if isfield(Hm0,'optstr')
    filename = strcat(filename,Hm0.optstr);
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% initialize modedata structure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % make copy of Hm0 m_structure
modedata        = Hm0;
% overwrite relevant properties
modedata.item   = shorttermdist(end);
modedata.label  = shorttermdist(end);
if strcmp(modedata.item,'H')
    modedata.bins = 2.0*Hm0.bins;
elseif strcmp(modedata.item,'C')
    modedata.bins = 1.5*Hm0.bins;
elseif strcmp(modedata.item,'T')
    modedata.bins = 1.1*Hm0.bins;
else
    error('Unknown item')
end

modedata.time   = [];
modedata.data   = [];
modedata.rho    = rho;
modedata.slope  = slope;
modedata.yrs    = yrs;
modedata.itemInfo = itemInfo;

% additional properties specific for mode data
modedata.ShortTermDist = shorttermdist;
modedata.Parameter     = parm;
modedata.Nyrs          = Nyrs;

if exist('WL','var')
    % time-varying water level given in input..
    modedata.WaterDepth = abs(WL.xyz(3)) + mean(WL.data(:,1),'omitnan');
    if strcmp(parm,'C') || strcmp(parm,'T')
        modedata.vref      = WL.vref;
        if ~strcmp(WL.vref,Hm0.vref)
            warning('WL.vref is %s while specified depth (Hm0.xyz(3)) is relative to %s.',WL.vref,Hm0.vref);
        elseif WL.xyz(3)~=Hm0.xyz(3)
            error('WL.vref (%s) equals Hm0.vref (%s), yet WL.xyz(3)=%.2f is different from Hm0.xyz(3)=%.2f',WL.vref,Hm0.vref,WL.xyz(3),Hm0.xyz(3));
        end
    else
        modedata.vref = '';
    end
else
    % constant water depth taken as Hm0.xyz(3)
    modedata.WaterDepth = abs(Hm0.xyz(3));
    if ~strcmpi(Hm0.vref,'MSL')
        warning('The water depth used as input to short-term distributions is relative to %s, not MSL. This may affect distributions in shallow water',Hm0.vref)
    end
    if strcmp(parm,'C') || strcmp(parm,'T')
        modedata.vref      = 'SWL';
    else
        modedata.vref = '';
    end
end
%
% vector for computing the CDF
X=(0:dx:ceil(max(Hm0.data)*3))';
DistIn.X = X;
%
% Interpolate to fill gaps (in omni-data only)
if interpol_hrs>0
    max_gap = round(interpol_hrs*3600/dur);
    inonan  = find(~isnan(Hm0.data(:,1)));
    tmp     = inonan(diff(inonan)>1 & diff(inonan)<=max_gap)+1;
    tmp     = meshgrid(tmp,1:max_gap)'+meshgrid(0:max_gap-1,1:length(tmp));
    iint    = intersect(find(isnan(Hm0.data(:,1))),reshape(tmp,[],1));
    Hm0.data(iint,1) = interp1(Hm0.time(inonan),Hm0.data(inonan,1),Hm0.time(iint),interpol_mthd);
    % Assume associated wave periods are missing for the same indices and
    % interpolate wave periods
    T02.data(iint,1) = interp1(T02.time(inonan),T02.data(inonan,1),T02.time(iint),interpol_mthd);
    for i=1:length(itemNames)
        if ~strcmpi(itemNames{i},{'Hm0','T02'})
            if itemDir(i)==0
                str = sprintf('%s.data(iint,1) = interp1(%s.time(inonan),%s.data(inonan,1),%s.time(iint),''%s'');',...
                    itemNames{i},itemNames{i},itemNames{i},itemNames{i},interpol_mthd);
            elseif itemDir(i)==1 % force nearest interpolation for directional items
                str = sprintf('%s.data(iint,1) = interp1(%s.time(inonan),%s.data(inonan,1),%s.time(iint),''nearest'');',...
                    itemNames{i},itemNames{i},itemNames{i},itemNames{i});
            end
            eval(str);
        end
    end
else
    iint = [];
end
%
% Decluster time series into independent events if iPeaks (storm indices)
% have not been provided as input to m_mode
if ~exist('iPeaks','var')
    eval(['declstrItem=',decluster_item,';']);
    if use_PDS
        iPeaks = find(PartialDurationSeries(declstrItem.time,declstrItem.data(:,1),intereventtime/24,intereventlevel));
        % define additional columns in iPeaks for compatibility
        iPeaks(:,2:5) = NaN;
        iPeaks(1,2) = 1;
        for r=1:size(iPeaks,1)-1      % loop over rows (all omni events)
            [~,tmp] = min(declstrItem.data(iPeaks(r,1):iPeaks(r+1,1),1));
            iPeaks(r,3)   = tmp+iPeaks(r,1)-1;
            iPeaks(r+1,2) = iPeaks(r,3)+1;
        end
        iPeaks(end,3) = length(declstrItem.time);
        iPeaks(:,4:5) = iPeaks(:,2:3);
    else
        iPeaks = m_decluster(declstrItem.time,declstrItem.data(:,1),'intereventtime',intereventtime,...
            'intereventlevel',intereventlevel,'IncludeStormTimeSteps',1,'StormCutOff',StormCutOff);
    end
    clear declstrItem
end
%
if size(iPeaks,2)<6
    iPeaks(:,6) = 1; % all events are analysed
end
%
% set size of output
nr = size(iPeaks,1);    % number of omni/all-year events
nc = size(Hm0.data,2);  % number of directional/monthly subdivisions (incl. omni)
Modes(nr,nc).Time = []; % init Modes structure
%
% Define lower Hm0 limits (to speed up computation by skipping small
% storms)
LowerLimit = NaN(1,nc);
if isempty(LowerLimitPrct)
    for c=1:nc; LowerLimit(c) = 0; end
else
    tmp = Hm0.data; tmp(isnan(tmp)) = 0;
    tmp = sort(tmp,'descend');
    for c=1:nc
        LowerLimit(c) = tmp(find(tmp(:,c)<=LowerLimitPrct(c)*max(tmp(:,c)),1,'first'),c);
    end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop through events and set event index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
events(nr).time = [];
events(nr).peak = [];
for r=1:nr      % loop over rows (all omni events)
    events(r).time = Hm0.time(iPeaks(r,1));
    events(r).peak = Hm0.data(iPeaks(r,1));
    events(r).indx = iPeaks(r,4:5); % use duration defined by Storm-cutoff
end
modedata.time = [events.time]';
for r=1:nr; for c=1:nc; Modes(r,c).Mode = NaN; end; end  % Fill Modes.Mode with NaNs


if return_only_maximum_mode>0 && return_only_maximum_mode<nr
    [~,imx] = sort([events.peak],'descend');
    events = events(imx(1:return_only_maximum_mode)); nr=return_only_maximum_mode; clear Modes
end

eventcount = zeros(1,nc);
for r=find([events.peak]>min_Hm0_peak & iPeaks(:,6)'==1)
    for c=1:nc
        % Duration defined by storm-cutoff. See m_decluster documentation
        index = intersect((events(r).indx(1):events(r).indx(2))',find(~isnan(Hm0.data(:,c))));
        
        % % reduce to time steps defined by StormCutOff
        % x = Hm0.data(index,c)/max(Hm0.data(index,c));
        % [~,ipeak] = max(x);
        % i0=find(x(1:ipeak)  <StormCutOff(1),1,'last');          if isempty(i0); i0=1; end
        % i1=find(x(ipeak:end)<StormCutOff(2),1,'first')+ipeak-1; if isempty(i1); i1=length(x); end
        % index = index(i0:i1);
        
        if ~isempty(index) && max(Hm0.data(index,c))>LowerLimit(c)
            eventcount(c) = eventcount(c)+1;
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Compute mode
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Modes(r,c).Hm0_peak = max(Hm0.data(index,c));
            imx = index(Hm0.data(index,c)==max(Hm0.data(index,c)));
            if numel(imx)>1 % in case same value is reached twice
                ind = [imx'-1;imx';imx'+1];
                [~,imx2] = max(sum(Hm0.data(ind)));
                imx = imx(imx2); clear imx2
            end
            Modes(r,c).Time     = Hm0.time(imx);
            Modes(r,c).StartTime= Hm0.time(index(1));
            Modes(r,c).EndTime  = Hm0.time(index(end));
            DistIn.Time         = Hm0.time(index);
            DistIn.Hm0          = Hm0.data(index);
            DistIn.T02          = T02.data(index);
            DistIn.N            = dur./DistIn.T02;
            %
            for i=1:length(itemNames)
                eval(['DistIn.',itemNames{i},' = ',itemNames{i},'.data(index);']);
            end
            %
            if exist('WL','var')
                % surge added to water depth. Triggers crest elevation calculation
                DistIn.d = abs(WL.xyz(3)) + WL.data(index);
            else
                % Constant water depth. Triggers crest height calculation
                DistIn.d = ones(size(DistIn.Hm0))*modedata.WaterDepth;
            end
            %
            % Call short term distribution function to return CDF
            DistOut = m_short_term_distributions_vFLD(shorttermdist,DistIn,1);
            CDF = DistOut.CDF;
            if ~isreal(CDF),CDF=real(CDF);end %SJA added because of A&P trough values
            nQ = (1-CDF)*DistIn.N;      % for eqv storm fit
            Pmax    = exp(-nQ);
            pmax    = diff(Pmax)./dx;
            Xmp     = interp1q(Pmax,X,exp(1)^-1);
            Modes(r,c).Mode = Xmp;
            %
            % transfer from crest heights to crest elevations
            if exist('WL','var') && (strcmp(parm,'C') || strcmp(parm,'T'))
                nQ_vref = zeros(size(X));
                for j=1:length(DistIn.Hm0)
                    nQ_vref = nQ_vref + interp1q(X+DistIn.d(j)-abs(WL.xyz(3)),(1-CDF(:,j))*DistIn.N(j),X);
                end
                Xmp_vref = interp1q(exp(-nQ_vref),X,exp(1)^-1);
                vref_str = [' (',num2str(Xmp_vref,'%5.2f'),'m',WL.vref,')'];
                Modes(r,c).Mode_vref = Xmp_vref;
            else
                vref_str = '';
                Modes(r,c).Mode_vref = [];
            end
            %
            % calculate sea state weight factors...
            Pmax0 = exp(-((1-CDF).*meshgrid(DistIn.N',1:length(DistIn.X))));
            wfac  = NaN(size(index));
            for i=1:length(index)
                tmpXmp  = interp1q(prod(Pmax0(:,setdiff(1:length(index),i)),2),X,exp(1)^-1); % Xmp with seastate(index(i)) excluded
                wfac(i) = Xmp-tmpXmp;
            end
            T02wmean = nansum(DistIn.T02.*wfac/sum(wfac(~isnan(DistIn.T02))));
            
            %%% Fit Equivalent (to distribution relative to SWL) %%%%%%%%%%%%%%%%%%
            if ~isempty(StartGuess)
                start_point = StartGuess(r,:);
            else
                start_point = [Modes(r,c).Hm0_peak log(max([(DistIn.Time(end)-DistIn.Time(1))*24*3600,dur])/mean(DistIn.T02)/2)];
            end
            EqvDistInput = DistIn;
            if strcmp(shorttermdist,'Forristall_C')
                EqvDistInput.alfa = DistOut.alfa_peak;
                EqvDistInput.beta = DistOut.beta_peak;
            end
            if strcmp(shorttermdist,'Battjes_Groenendijk_H')
                EqvDistInput.Htr = DistOut.Htr;
            end
            if exist('WL','var')
                EqvDistInput.d  = abs(WL.xyz(3))+nansum(DistIn.WL.*wfac/sum(wfac(~isnan(DistIn.WL)))); % use the weighted average mean WL for equivalent storm
            else
                EqvDistInput.d  = abs(Hm0.xyz(3));
            end
            
            summary_str = sprintf('Storm no: %i,%i - %s. Peak Hm0=%5.2f. %smp=%5.2fm%s.',...
                r,c,datestr(Modes(r,c).Time(1),'dd-mm-yyyy'),max(DistIn.Hm0),parm,Xmp,vref_str);
            %
            pmax_eqvBox=[];  eqGauss=[];
            if ismember(1,output_types)
                % fit equivalent box storm
                [eqBox, model] = FitEqvStorm(X,pmax,start_point,shorttermdist,EqvDistInput);
                [~, pmax_eqvBox] = model(eqBox);
                eqBox(2) = round(exp(eqBox(2)));
                summary_str = [summary_str,sprintf(' Hm0_eqv=%5.2fm, N_eq=%.0f',eqBox)];
                %Xmp_eq=X(pmax_eq==max(pmax_eq))+dx/2;
                %
                % fit Gumbel distribution to Pmax to get N
                [eqGumb, model] = FitGumbel(X,Pmax,start_point,shorttermdist,EqvDistInput);
                [~, Pmax_gb] = model(eqGumb); pmax_gb = diff(Pmax_gb)/dx;
                eqGumb(2) = round(exp(eqGumb(2)));
                %
                Modes(r,c).Hm0_eqv      = eqBox(1);
                Modes(r,c).N_eqv        = eqBox(2);
                Modes(r,c).eqGumb_Xmp   = eqGumb(1);
                Modes(r,c).eqGumb_N     = eqGumb(2);
                
            end
            if ismember(2,output_types)
                % fit equivalent Gauss storm (to be used with J-EVA)
                [eqGauss, model] = FitSTDStorm(X,Pmax,[start_point(1),start_point(2)],shorttermdist,EqvDistInput,mnsz_lnlnP);
                [~,pmax_Gauss] = model(eqGauss);
                Modes(r,c).Hm0peq       = eqGauss(1);
                Modes(r,c).lnsigmaeq    = eqGauss(2);
                eqGauss(2) = exp(eqGauss(2));
                summary_str = [summary_str,sprintf(' Hm0peq: %5.2fm, sigma_eq: %.0f',eqGauss)];
            end
            fprintf('%s\n',summary_str);
            
            Modes(r,c).TSvalues = rmfield(DistIn,'X');
            Modes(r,c).interpol_frac= length(intersect(iint,events(r).indx(1):events(r).indx(2)))/(diff(events(r).indx)+1);
            %
            % Weighted average of associated parameters
            Modes(r,c).T02wmean = T02wmean; % weighted mean T02
            for i=1:length(itemNames)
                if ~strcmpi(itemNames{i},{'Hm0','T02'})
                    if eval(['any(~isnan(DistIn.',itemNames{i},'))'])
                        if itemDir(i)==0
                            eval(['Modes(r,c).',itemNames{i},'wmean = nansum(DistIn.',itemNames{i},...
                                '.*wfac/sum(wfac(~isnan(DistIn.',itemNames{i},'))));']);
                        elseif itemDir(i)==1
                            eval(['Uwmean = nansum(DistIn.Hm0.*sin(DistIn.',itemNames{i},'.*pi/180).*wfac/sum(wfac(~isnan(DistIn.',itemNames{i},'))));']);
                            eval(['Vwmean = nansum(DistIn.Hm0.*cos(DistIn.',itemNames{i},'.*pi/180).*wfac/sum(wfac(~isnan(DistIn.',itemNames{i},'))));']);
                            dir = atan2(Uwmean,Vwmean).*180/pi;if dir<0,dir=dir+360;end
                            eval(['Modes(r,c).',itemNames{i},'wmean = dir;']);
                        end
                    end
                end
            end
            
            if strcmp(shorttermdist,'Forristall_C')
                Modes(r,c).alfa_peak = DistOut.alfa_peak;
                Modes(r,c).beta_peak = DistOut.beta_peak;
            end
            if strcmp(shorttermdist,'Glukhovskiy_H')
                Modes(r,c).kappa_peak = DistOut.kappa_peak;
            end
            if strcmp(shorttermdist,'Weibull_H')
                modedata.alpha = alpha;
            end
            
            if plotflag==1
                index0 = (max(1,(index(1)-round(dur/300))):min(size(Hm0.data,1),(index(end)+round(dur/300))))'; % index only for plots (storm +-1/2day)
                %%% Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                close
                clear lh str1
                subplot(6,2,1:4),hold on, grid on
                bar(Hm0.time(index0),Hm0.data(index0,c),'EdgeColor',Hm0.ColorOrder(3,:),'FaceColor','w')
                lhdl=bar(DistIn.Time,DistIn.Hm0,'EdgeColor',Hm0.ColorOrder(3,:),'facecolor',Hm0.ColorOrder(3,:));
                lstr = {'H_{m0}  [m]'};
                if ~isempty(iint)
                    bar(Hm0.time(iint),Hm0.data(iint,c),'EdgeColor',Hm0.ColorOrder(3,:),'facecolor',0.7*ones(1,3))
                end
                lhdl(2)=plot(DistIn.Time,DistIn.T02,'^','color',Hm0.ColorOrder(2,:),'markersize',5,'markerfacecolor',Hm0.ColorOrder(2,:));
                lstr{2} = sprintf('T_{02} (w_{mean}=%.1fs)',Modes(r,c).T02wmean);
                plot(DistIn.Time([1,end]),T02wmean*[1,1],'-','color',Hm0.ColorOrder(2,:),'linewidth',1)
                tstr = sprintf('%s: H_{m0,peak}=%.2fm, %s_{mp}=%.2fm',...
                    datestr(Modes(r,c).Time,'dd-mmm-yyyy'),Modes(r,c).Hm0_peak,parm,Modes(r,c).Mode);
                
                if ~isempty(eqGauss)
                    % plot Gauss bell shape above 70% of peak value
                    x = eqGauss(2)*sqrt(-2*log(0.7))*linspace(-1,1,51);
                    y = eqGauss(1)*exp(-x.^2/2/eqGauss(2)^2);% pdfnorm(x,0,eqGauss(2)^2)/pdfnorm(0,0,eqGauss(2)^2);
                    plot(x*T02wmean/3600/24+sum(DistIn.Time.*wfac/sum(wfac)),y,'-k','LineWidth',1.2)
                    Modes(r,c).bellTime = x([1,end])*T02wmean/3600/24+sum(DistIn.Time.*wfac/sum(wfac));
                    tstr = [tstr,sprintf(', H_{m0,p,eq}=%.2fm, \\sigma_{eqN}=%.0f',eqGauss)];
                else
                    tstr = [tstr,sprintf(', H_{m0,eqv box}=%.2fm, N_{eqv}=%.0f',eqBox)];
                end
                
                if exist('T01','var')
                    lhdl(end+1)=plot(DistIn.Time,DistIn.T01,'d','color',Hm0.ColorOrder(5,:),'markersize',5,'markerfacecolor',Hm0.ColorOrder(5,:));
                    lstr{end+1} = sprintf('T_{01} (w_{mean}=%.1fs)',Modes(r,c).T01wmean);
                    plot(DistIn.Time([1,end]),Modes(r,c).T01wmean*[1,1],'-','color',Hm0.ColorOrder(5,:),'LineWidth',1);
                end
                if exist('Tp' ,'var')
                    lhdl(end+1)=plot(DistIn.Time,DistIn.Tp ,'v','color',Hm0.ColorOrder(4,:),'markersize',5,'markerfacecolor',Hm0.ColorOrder(4,:));
                    lstr{end+1} = sprintf('T_{p} (w_{mean}=%.1fs)',Modes(r,c).Tpwmean);
                    plot(DistIn.Time([1,end]),Modes(r,c).Tpwmean*[1,1],'-','color',Hm0.ColorOrder(4,:),'LineWidth',1);
                end
                legend(lhdl,lstr,'Location','SouthEast')
                xtick = floor(DistIn.Time(1)*2)/2:12/24:ceil(DistIn.Time(end)*2)/2;
                set(gca,'XTick',xtick,'XTickLabel',char(datestr(xtick(1),'dd-mm-yy HH:MM'),datestr(xtick(2:end),'HH:MM')),...
                    'FontSize',fontsize,'XLim',xtick([1,end])+[-0.1,0.1])
                %m_xticklabel_rotate([],45)
                % text(Modes(r,c).Time,Modes(r,c).Hm0_peak+1,['H_{m0,peak} = ',...
                %     num2str(Modes(r,c).Hm0_peak,'%5.2f'),'m. ',parm,'_{mp} = ',...
                %     num2str(Modes(r,c).Mode,'%5.2f'),'m.'],'FontSize',fontsize)
                % title(sprintf('%s: H_{m0,peak}=%.2fm, %s_{mp}=%.2fm, H_{m0,p,eq}=%.2fm, \\sigma_{eqN}=%.0f',...
                %     datestr(Modes(r,c).Time,'dd-mmm-yyyy'),Modes(r,c).Hm0_peak,parm,Modes(r,c).Mode,eqGauss(1),eqGauss(2)))
                title(tstr)
                
                for sb=1:3
                    cc=0;lhdl = []; lstr=[]; ilh=1; mk = {'^','d','v'}; cl=Hm0.ColorOrder([2,4,5],:);
                    switch sb
                        case 1; sbv=5:6;  items = {'WS','DSD'};       ndec=[1,1];   ylm=[];
                        case 2; sbv=7:8;  items = {'WL','CS'};        ndec=[2,2];   ylm=[];
                        case 3; sbv=9:10; items = {'MWD','PWD','WD'}; ndec=[0,0,0]; ylm=[0,360];
                    end
                    for it = 1:length(items)
                        if exist(items{it},'var')
                            cc=cc+1;
                            eval(['item = ',items{it},';']);
                            eval(['ts = [DistIn.Time,DistIn.',items{it},'];']);
                            eval(['wmean = Modes(r,c).',items{it},'wmean;']);
                            subplot(6,2,sbv),hold on,grid on
                            lhdl(cc) = plot(ts(:,1),ts(:,2),mk{it},'color',cl(it,:),'markersize',5,'markerfacecolor',cl(it,:));
                            lstr{cc} = sprintf('%s (w_{mean}=%s%s)',item.label,num2str(wmean,['%.',int2str(ndec(it)),'f']),item.unit);
                            plot(ts(:,1),ones(size(ts(:,1)))*wmean,'-','color',cl(it,:),'LineWidth',1);
                        end
                    end
                    if ~isempty(lhdl)
                        set(gca,'XTick',xtick,'XTickLabel',char(datestr(xtick(1),'dd-mm-yy HH:MM'),datestr(xtick(2:end),'HH:MM')),...
                            'FontSize',fontsize,'XLim',xtick([1,end]))
                        legend(lhdl,lstr,'Location','SouthEast');
                        if ~isempty(ylm); set(gca,'ytick',0:45:360,'ylim',ylm); end
                    end
                end
                
                % pdf comparison plot
                subplot(6,2,11)
                lh = plot(X(2:end)-dx/2,pmax   ,'color',Hm0.ColorOrder(3,:),'Marker','.');hold on
                lstr = {'Actual Storm'};
                if ~isempty(eqBox)
                    lh(end+1)=plot(X(2:end)-dx/2,pmax_eqvBox,'color',Hm0.ColorOrder(2,:),'LineWidth',1.2); lstr(end+1)={'Eqv. Box'};
                    lh(end+1)=plot(X(2:end)-dx/2,pmax_gb,'color',Hm0.ColorOrder(4,:)); lstr(end+1)={'Gumbel'};
                end
                if ~isempty(eqGauss)
                    lh(end+1)=plot(X(2:end)-dx/2,pmax_Gauss,'color',Hm0.ColorOrder(5,:)); lstr(end+1)={'Eqv. Gauss'};
                end
                ylm=get(gca,'YLim');ylm(1)=0;
                set(gca,'FontSize',fontsize,'YLim',ylm)
                if ~isempty(eqGauss)
                    text(0.01*max(X),0.30*ylm(2),['H_{m0,p,eq}=',num2str(eqGauss(1),'%5.2f'),'m'],'FontSize',0.6*fontsize,'BackgroundColor','w')
                    text(0.01*max(X),0.15*ylm(2),['\sigma_{eq}=',int2str(eqGauss(2))],'FontSize',0.6*fontsize,'BackgroundColor','w')
                else
                    text(0.01*max(X),0.60*ylm(2),['H_{m0,eqv box}=',num2str(eqBox(1),'%5.2f'),'m'],'FontSize',0.6*fontsize,'BackgroundColor','w')
                    text(0.01*max(X),0.45*ylm(2),['N_{eqv}=',int2str(eqBox(2))],'FontSize',0.6*fontsize,'BackgroundColor','w')
                    text(0.01*max(X),0.30*ylm(2),[parm,'_{mp,Gumbel}=',num2str(eqGumb(1),'%5.2f'),'m'],'FontSize',0.6*fontsize,'BackgroundColor','w')
                    text(0.01*max(X),0.15*ylm(2),['N_{Gumbel}=',int2str(eqGumb(2))],'FontSize',0.6*fontsize,'BackgroundColor','w')
                end
                ylabel(['p(',parm,'_{max}|storm)']),%legend('Actual Storm','Eq. Storm')
                xlabel([parm,'_{max} [m]'])
                axis tight
                legend(lh,lstr)
                
                % ln(-ln(P)) comparison plot
                subplot(6,2,12),hold on,grid on, box on
                set(gca,'FontSize',fontsize)
                plot(X,-log(-log(Pmax))   ,'color',Hm0.ColorOrder(3,:),'Marker','.'),hold on
                if ~isempty(eqGauss)
                    plot(X,-log(-log([0;cumsum(pmax_Gauss)/sum(pmax_Gauss)])),'color',Hm0.ColorOrder(5,:))
                end
                axis tight
                xlabel([parm,'_{max} [m]']),ylabel('-ln(-ln(P))')
                
                if ~isempty(OutputFolder)
                    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 32 42])
                    print('-dpng','-r150',strcat(OutputFolder,filename,'_',datestr(Modes(r,c).Time,'yyyy-mm-dd_HH'),Hm0.figure)), close
                end
                if exist('shortpngfilename','var')
                    print('-dpng','-r150',shortpngfilename)
                end
                %%%%
                
            end
            
            clear nQ_vref
        else
            Modes(r,c).Mode      = NaN;
            Modes(r,c).Mode_vref = NaN;
        end % if index vector is not empty
    end     % loop over columns of Events (=no. of directions+1)
end         % loop over rows of Events    (=nyrs for AMP extraction)
%
% Load peak time and modes into m_structure format
if strcmp(modedata.vref,Hm0.vref) % HFH renamed from '_MSL' to more general '_vref'
    modedata.data = reshape([Modes.Mode_vref],size(Modes));
else
    modedata.data = reshape([Modes.Mode],size(Modes));
end
%
% Also save additional mode data information in Modes sub-structure
modedata.Modes = Modes;
modedata.filename = filename;
modedata.declustered_by = decluster_item;
%
if remove_TS % remove TSvalues field in order to reduce file size
    modedata.Modes=rmfield(modedata.Modes,'TSvalues');
end
%
% Save results if requested
if ~isempty(OutputFolder)
    save(strcat(OutputFolder,filename,'.mat'),'modedata');
end
%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [estimates, model] = FitEqvStorm(X, p, start_point, ShortTermDist, DistIn)
% Fit equivalent box-shaped storm
options = optimset;
options.MaxIter=1e4;
options.MaxFunEvals=1e4;
options.ShortTermDist = ShortTermDist;
model = @Eq;
estimates = fminsearch(model, start_point, options);

    function [sse, pfit] = Eq(params)
        Hm0  = params(1);
        N    = exp(params(2));
        DistIn.X   = X;
        DistIn.Hm0 = Hm0;
        
        DistOut     = m_short_term_distributions_vFLD(options.ShortTermDist,DistIn,1);
        P           = DistOut.CDF.^N;
        pfit        = diff(P)/(X(2)-X(1));
        ErrorVector = pfit - p;
        sse         = sum(ErrorVector .^ 2);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [estimates, model] = FitGumbel(X, P, start_point, ShortTermDist, DistIn)
% Fit Gumbel distribution to Pmax
options = optimset;
options.MaxIter=1e4;
options.MaxFunEvals=1e4;
model = @Gumb;
DistIn.Xmx = X;
estimates = fminsearch(model, start_point, options);

    function [sse, Pfit] = Gumb(params)
        DistIn.X = params(1);
        DistIn.N = exp(params(2));
        
        Pfit = m_short_term_distributions_vFLD(ShortTermDist,DistIn,3);
        ErrorVector = Pfit - P;
        sse         = sum(ErrorVector .^ 2);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [estimates, model] = FitSTDStorm(X, P, start_point, ShortTermDist, DistIn, mnsz_lnlnP)
% Fit peak and standard deviation of Gauss bell shaped storm
options = optimset;
options.MaxFunEvals=100000;
options.ShortTermDist = ShortTermDist;
model = @Fit_STD;
d = DistIn.d; % d is scalar
p = diff(P)./diff(X);
estimates = fminsearch(model, start_point, options);
    function [sse, pfit, Pfit] = Fit_STD(params)
        Hm0p        = params(1);
        STD         = exp(params(2));
        npts        = 51;
        N           = 2*sqrt(-2*STD^2*log(0.7));                     % Truncate Gauss bell at 70% of peak
        Hm0         = Hm0p*exp(-linspace(-N/2,N/2,npts).^2/2/STD^2); % Gauss bell shape
        DistIn.Hm0  = Hm0;
        DistIn.N    = ones(npts,1)*N/npts;
        DistIn.d    = d*ones(size(Hm0));                             % create a vector of identical depths
        DistOut     = m_short_term_distributions_vFLD(options.ShortTermDist,DistIn,1);
        nQ          = (1-DistOut.CDF)*DistIn.N;      % for eqv storm fit
        Pfit        = exp(-nQ);
        pfit        = diff(Pfit)./diff(X);
        if mnsz_lnlnP
            ii = Pfit>eps & Pfit<1-eps & P>eps & P<1-eps;
            ErrorVector = log(-log(Pfit(ii))) - log(-log(P(ii)));
        else
            ErrorVector = pfit - p;
        end
        sse         = sum(ErrorVector .^ 2);
    end
end
%
%
function assign_argin(argum)
for i=1:2:length(argum)
    assignin('caller',argum{i},argum{i+1});
end
end
