function X = m_resample(X,varargin)
% modifed by jwo
% Resample X to an equidistant time series with time step of dt (ttt(4))
% using nearest value in time and replacing data values longer apart than
% maxgap (ttt(4)/2) with nan. E.g. '10min values @ 1h intervals' or '10min
% values @ 1min intervals'.
%
% SYNTAX:    X = m_resample(X)
%            X = m_resample(X,'ttt',ttt)
%            m_resample(X,'maxgap',60)
%            X = m_resample(X,'dt_avg',120) 
%
% INPUT:     X          = m_structure
%            ttt        = ttt ['begin date' 'end date' 'averging time (min)']
%            dt_avg     = Averging time [min] (default = X.ttt(3))
%            dt_int     = time interval [min] (default = X.ttt(4))
%            maxgap     = Max gap in time [min] (default = X.ttt(4)/2)
%            vector     = Vector of magnitude for directional resampling
%            gap_fill   = gap filling using interp1 [default = 0]
%
% OUTPUT:    Equidistant time series with dt time step
%
% Copyright DHI ©, Denmark, All rights reserved
% 2016-02-02 - Patrick Dich Grode, pdg@dhigroup.com
% 2018-04-03 - PDG: added gap_fill

% Dev to consider (from HEWR):
% TT = timetable(Obs.ED1f.datetime,Obs.ED1f.data);
% TT = smoothdata(TT,'movmean',hours(3));
% Obs.ED1f.data = TT.Var1;

% Varargins
dt_avg = [];
dt_int = [];
maxgap = [];
Y      = [];
gap_fill = 0;
if ~isempty(varargin)
    i = 1;
    while i <= length(varargin)
        if     strcmp(varargin{i},'ttt')
            X.ttt = varargin{i+1};                                 i = i+1;
        elseif strcmp(varargin{i},'dt_avg')
            dt_avg = varargin{i+1};                                i = i+1;
        elseif strcmp(varargin{i},'dt_int')
            dt_int = varargin{i+1};                                i = i+1;
        elseif strcmp(varargin{i},'maxgap')
            maxgap = varargin{i+1};                                i = i+1;
        elseif strcmp(varargin{i},'vector')
            Y = varargin{i+1};                                     i = i+1;
        elseif strcmp(varargin{i},'gap_fill')
            gap_fill = 1;                                          i = i+0;
        else
            error('Unknown argument: %s',num2str(varargin{i}))
        end
        i = i+1;
    end
end

% Averaging time [min]
if isempty(dt_avg), dt_avg = X.ttt(3); end

% Time interval [min]
if isempty(dt_int), dt_int = X.ttt(4); end

% Maxgap in time [min]
if isempty(maxgap), maxgap = dt_int/2; end

% if isstruct(Y) % PDG: Why?
%     Y = Y.data;
% end

% Round t_min and t_max
t_min = X.ttt(1)-rem(X.ttt(1),dt_int/24/60); % floor
t_max = X.ttt(2)-rem(X.ttt(2),dt_int/24/60); if rem(X.ttt(2),dt_int/24/60) ~= 0, t_max = t_max + dt_int/24/60; end % ceil

% Check time interval
dt_mode = mode(diff(X.time));
if dt_int < 0.99*dt_mode
    warning('Given time interval is less than 99% of most frequent time interval (i.e. the resampled time series contain many nans)')
end

% Remove nans (if nearest point is nan and 2. is not)
ndatadims = ndims(X.data)-1;
if ndatadims == 1 % old for ndims==1
    Ifin = isfinite(X.data);
    X.time = X.time(Ifin);
    X.data = X.data(Ifin);
    clear Ifin
end % skip this step if multiple dimensions

% Create equidistant time vector ranging ttt
time = [t_min:dt_int/60/24:t_max]';

% Find indices of nearest times in sample vector
% Use extrap to edges & to avoid nans in output
if isempty(X.time) || length(X.time) < 2
    disp('m_resample: X.time was empty!'), return
else
    It = interp1(X.time,1:length(X.time),time,'nearest','extrap');
end

% select data
if ndatadims == 1
    % Create data vector associated with time vector
    data = X.data(It);
    
    % Reject data >= than single(maxgap) apart in time
    data(single(abs(time-X.time(It))) >= single(maxgap/60/24)) = NaN; clear It
else
    data = X.data(It,:);
    m = single(abs(time-X.time(It))) >= single(maxgap/60/24);
    data(m,:) = NaN;
    if ndatadims>2
        s = size(X.data);
        snew = s; snew(1) = sum(m);
        data = reshape(data,snew);
    end
end


% Output
clear X.time X.data;  X.time = time; X.data = data; clear data
if ismember('datetime',fieldnames(X)) % HEWR: add passing through of .datetime
    clear X.datetime; X.datetime = datetime(X.time,'ConvertFrom','datenum');
end

% Running average of data if new averaging time is given
dt_avg=60;
try dt_avg_X = X.ttt(3); catch, dt_avg_X = nan; end
if dt_avg ~= dt_avg_X && ~isnan(dt_avg)
    disp('Time averaging')
    % Vector average if dir
    if X.isdir
        
        % Y (possibly not same time vector as X!)
        if ~isempty(Y)
            It = interp1(Y.time,1:length(Y.time),time,'nearest','extrap');
            data = Y.data(It);
            data(single(abs(time-Y.time(It))) >= single(maxgap/60/24)) = NaN; clear It
            clear Y.time Y.data;  Y.time = time; Y.data = data; clear time data
        end

        u = X;u.isdir = 0;
        v = X;v.isdir = 0;
        if ~isempty(Y)
            [u.data,v.data] = spddir2uv(Y.data,X.data);  u.isdir = 0;  v.isdir = 0;
            u     = m_resample(u,'dt_avg',dt_avg);
            v     = m_resample(v,'dt_avg',dt_avg);
            [~,X.data] = uv2spddir(u.data,v.data); X.isdir = 1; % return
        else
            warning('Directions should be vector averaged')
            [u.data,v.data] = spddir2uv(ones(size(X.data)),X.data);  u.isdir = 0;  v.isdir = 0;
            u     = m_resample(u,'dt_avg',dt_avg);
            v     = m_resample(v,'dt_avg',dt_avg);
            [~,X.data] = uv2spddir(u.data,v.data); X.isdir = 1; % return
        end
        
    else
        
        % Check averaging time
        if rem(dt_avg,dt_int) > 0
            error('Averaging time (ttt(3)) must be a multiple of time interval (ttt(4))')
        end
        
        % Check data coverage
        if dt_avg < 0.50*dt_mode
            warning('Data coverage for given averaging time is less than 50% of most frequent time interval')
        end
        
        % Running average
        N = length(X.time);
        di1 = floor(dt_avg/dt_int/2);
        di2 = ceil(dt_avg/dt_int/2);
        
        % display('Running averaging requesting >= 50% finite values')
        M = nan(N,di1+di2); % PDG: removed +1 in col size
        for iM = 1:di1+di2
            if iM <= di1
                M(1:end+iM-di1-1,iM) = X.data(di1-iM+2:end);
            elseif iM == di1+1
                M(:,iM) = X.data;
            else
                M(iM-di1:end,iM) = X.data(1:end-iM+di1+1);
            end
        end
        
        % Gap filling
        if gap_fill
            for j = 1:size(M,2)
                nanx = isnan(M(:,j)); t = 1:numel(M(:,j)); M(nanx,j) = interp1(t(~nanx), M(~nanx,j), t(nanx));
            end
        end
        
        % M(isnan(M(:,di1+1)),:) = NaN; % PDG: commented out, central value defines nan, why?
        data = nanmean(M,2);
        data(sum(isfinite(M')) < di2) = NaN; % Request >= 50% finite values
        
        clear X.data;  X.data = data; clear data
    end

else
    disp('No time averaging performed.')
end

% Update heading with resampled ttt - dt_avg
X.ttt = [t_min t_max dt_avg dt_int];
X = m_heading(X);

% Updated datetime if existed before
if ismember('datetime',fieldnames(X))
    X.datetime = datetime(X.time,'ConvertFrom','datenum');
end

end
