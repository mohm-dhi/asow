function h = m_smoothhist2D_mohm(X,lambda,nbins,outliercutoff,plottype,cmap,bin2,unit,unitY,edges1,edges2,densityPower, marker_size)
% SMOOTHHIST2D Plot a smoothed histogram of bivariate data.
%   SMOOTHHIST2D(X,LAMBDA,NBINS) plots a smoothed histogram of the bivariate
%   data in the N-by-2 matrix X.  Rows of X correspond to observations.  The
%   first column of X corresponds to the horizontal axis of the figure, the
%   second to the vertical. LAMBDA is a positive scalar smoothing parameter;
%   higher values lead to more smoothing, values close to zero lead to a plot
%   that is essentially just the raw data.  NBINS is a two-element vector
%   that determines the number of histogram bins in the horizontal and
%   vertical directions.
%
%   SMOOTHHIST2D(X,LAMBDA,NBINS,CUTOFF) plots outliers in the data as points
%   overlaid on the smoothed histogram.  Outliers are defined as points in
%   regions where the smoothed density is less than (100*CUTOFF)% of the
%   maximum density.
%
%   SMOOTHHIST2D(X,LAMBDA,NBINS,[],'surf') plots a smoothed histogram as a
%   surface plot.  SMOOTHHIST2D ignores the CUTOFF input in this case, and
%   the surface plot does not include outliers.
%
%   SMOOTHHIST2D(X,LAMBDA,NBINS,CUTOFF,'image') plots the histogram as an
%   image plot, the default.
%
%   Example:
%       X = [mvnrnd([0 5], [3 0; 0 3], 2000);
%            mvnrnd([0 8], [1 0; 0 5], 2000);
%            mvnrnd([3 5], [5 0; 0 1], 2000)];
%       smoothhist2D(X,5,[100, 100],.05);
%       smoothhist2D(X,5,[100, 100],[],'surf');
%
%   Reference:
%      Eilers, P.H.C. and Goeman, J.J (2004) "Enhancing scaterplots with
%      smoothed densities", Bioinformatics 20(5):623-628.

%   Copyright 2009 The MathWorks, Inc.
%   Revision: 1.0  Date: 2006/12/12
%
%   Requires MATLAB® R14.

if nargin < 4 || isempty(outliercutoff), outliercutoff = .05; end
if nargin < 5, plottype = 'image'; end
if nargin < 9, unitY = []; end %unit Y assumed = unit X
if nargin < 12, densityPower = 0.5; end
                

% if only one bin size is given it will be used on both axis.
if length(bin2) == 1
    bin2(1:2) = bin2;
end

if exist('edges1','var')
    ctrs1 = edges1(1:end-1) + .5*diff(edges1);
    ctrs2 = edges2(1:end-1) + .5*diff(edges2);
else
    minx = min(X(:,1));
    maxx = max(X(:,1));
    minx(2) = min(X(:,2));
    maxx(2) = max(X(:,2));
    
    edges1 = linspace(minx(1), maxx(1), nbins(1)+1);
    ctrs1 = edges1(1:end-1) + .5*diff(edges1);
    edges1 = [-Inf edges1(2:end-1) Inf];
    edges2 = linspace(minx(2), maxx(2), nbins(2)+1); %pje
    ctrs2 = edges2(1:end-1) + .5*diff(edges2);
    edges2 = [-Inf edges2(2:end-1) Inf];
end

% F = m_histcn([X(:,1) X(:,2)],  edges1, edges2);
% % if there are data in X that is equal edges(end), then the size of F
% % increased by 1, this is ten removed again below.
% if size(F,1)>(length(edges1)-1)
%     F = F(1:end-1,:);
% end
% if size(F,2)>(length(edges2)-1)
%     F = F(:,1:end-1);
% end
% PDG: update adopted from aux00_smoothhist2D (jeva - hfh):
[Fcount,xedges,yedges,binX,binY] = histcounts2(X(:,1),X(:,2),edges1,edges2); % histcounts2 includes data that fall on the edges..
F = Fcount;
% linind = sub2ind(size(F),binX,binY);

% PDG: This ruins/changes density and contours in smoe cases!
% F(F<1) = 1;
% F(F>m_quantile(F(:),0.9995))= m_quantile(F(:),0.9995);
F = F.^densityPower;
F = smooth1D(F,lambda);



% 
% [n,p] = size(X);
% bin = zeros(n,2);
% % Reverse the columns of H to put the first column of X along the
% % horizontal axis, the second along the vertical.
% [dum,bin(:,2)] = histc(X(:,1),edges1);
% [dum,bin(:,1)] = histc(X(:,2),edges2);
% bin(bin==0) = 1; %hack by pje
% H = accumarray(bin,1,round(nbins([2 1]))) ./ n;
% 
% % Eiler's 1D smooth, twice
% G = smooth1D(H,lambda);
% F = smooth1D(G',lambda)';
% % % An alternative, using filter2.  However, lambda means totally different
% % % things in this case: for smooth1D, it is a smoothness penalty parameter,
% % % while for filter2D, it is a window halfwidth
% % F = filter2D(H,lambda);
% F = F.*size(X,1);
% relF = F./max(max(F(:)));
% if outliercutoff > 0
%     outliers = (relF(nbins(2)*(bin(:,2)-1)+bin(:,1)) < outliercutoff);
% end

% nc = 256;
% colormap(hot(nc));
switch plottype
    case {'scatter'}
        nc = length(cmap);
        colormap(cmap)
    case {'image'}
        nc = length(cmap);
        colormap([1 1 1; cmap])% add white for background
        alpha = zeros(size(F));
        alpha(floor(nc.*relF) + 1 ~= 1) =1;
    otherwise
        colormap(cmap)
%         F(F<outliercutoff) = nan;
end

switch plottype
    case 'surf'
        h = surf(ctrs1,ctrs2,F./nanmax(nanmax(F))-1,'edgealpha',0);
        hold on
        view(2)
        shading flat
        %         if outliercutoff > 0
        %             plot(X(outliers,1),X(outliers,2),'.','MarkerEdgeColor',[.8 .8 .8]);
        %         end     
    case 'image'
        h = image(ctrs1,ctrs2,floor(nc.*relF) + 1,'alphadata',alpha);
        hold on
    case 'gridscatter'
        [x,y]=meshgrid(ctrs1,ctrs2);
%         F(F<=outliercutoff)=nan;
        h = scatter(x(:),y(:),1,F(:));
    case 'scatter'
         [x,y]=meshgrid(ctrs1,ctrs2);
%         C = interp2(x,y,F,X(:,1),X(:,2),'nearest',1); % ceil ,'nearest',1 - only if no smoothing!
        C = interp2(x,y,F',X(:,1),X(:,2),'spline',1);
        %CMaxLim = interp1(linspace(0,1,length(C(~isnan(C)))),sort(C(~isnan(C))),1); % PDG xq = 1 instead of 0.99?
        %C = min(C,CMaxLim);
         h = gsp(X(:,1),X(:,2),C,marker_size,cmap);
        % PDG: try log instead:
        %h = gsp(X(:,1),X(:,2),log(C),10,cmap);
end

% PDG quick and dirty:
t = colorbar('peer',gca);
if length(X) < 3000
    d = 1;
elseif length(X) < 30000
    d = 10;
elseif length(X) < 300000
    d = 100;
elseif length(X) < 3000000
    d = 1000;
else
    d = 10000;
end
switch plottype
    case {'image','scatter'}
        C = C.^(1/densityPower);
        d_scale = max(round((max(max(C,2)) - max(min(C),1))/d)*d/10,1);
        c_tick_label = unique([max(min(C),1) d_scale : d_scale :  max(max(C,2))]); % use 1 as min value and 2 as max value
        
% PDG: try log instead:
        
%         d_scale = max(round((max(max(log(C),2)) - max(min(log(C)),1))/log(d))*log(d)/10,1);  % c_tick_label = reshape([1 3]'*10.^(0:7),1,[]);
%         c_tick_label = round(exp(unique([max(min(log(C)),1) d_scale : d_scale :  max(max(log(C),2))]))); % use 1 as min value and 2 as max value
        
        while length(c_tick_label)>30
            d = 10^(log10(d)-1);
            d_scale = max(round((max(max(C,2)) - max(min(C),1))/d)*d/10,1);
            c_tick_label = unique([max(min(C),1) d_scale : d_scale :  max(max(C,2))]); % use 1 as min value and 2 as max value
        end
        matver = version('-release');
        if str2num(matver(1:4))<2014
            c_tick = unique(nc/(max(max(C,2))-max(min(C),1)).*(c_tick_label-max(min(C),1))+1); % use 1 as min value and 2 as max value
        else
            c_tick = linspace(0,1,length(c_tick_label));
        end
        set(t,'ytick',c_tick.^(densityPower))
        set(t,'yticklabel',c_tick_label) % scale from color scale to date points scale
    case 'surf'
        c_scale = (get(t,'ytick')+1).*nanmax(nanmax(F));
        dScale = round(mean(diff(c_scale))/d)*d;
        %c_scale = d*round(c_scale / d);
        c_scale = (d*floor(c_scale(1)/d)):dScale:(d*ceil(c_scale(end)/d));
        set(t,'ytick',c_scale./nanmax(nanmax(F))-1)
        set(t,'yticklabel',c_scale) % scale from color scale to date points scale
        %         set(t,'yticklabel',(get(t,'ytick')+1).*nanmax(nanmax(F)));
end
if isnumeric(unitY) && bin2(1)==bin2(2)
    set(get(t,'ylabel'),'String',['Number of data points in each ' num2str(bin2(1)) ' ' unit ' bin']);
else
    set(get(t,'ylabel'),'String',['Number of data points in each ' num2str(bin2(1)) ' ' unit ,' x ' num2str(bin2(2)) ' ' unitY ' cell']);
end



%-----------------------------------------------------------------------------
function Z = smooth1D(Y,lambda)
[m,n] = size(Y);
E = eye(m);
D1 = diff(E,1);
D2 = diff(D1,1);
P = lambda.^2 .* D2'*D2 + 2.*lambda .* D1'*D1;
Z = (E + P) \ Y;
% This is a better solution, but takes a bit longer for n and m large
% opts.RECT = true;
% D1 = [diff(E,1); zeros(1,n)];
% D2 = [diff(D1,1); zeros(1,n)];
% Z = linsolve([E; 2.*sqrt(lambda).*D1; lambda.*D2],[Y; zeros(2*m,n)],opts);


%-----------------------------------------------------------------------------
function Z = filter2D(Y,bw)
z = -1:(1/bw):1;
k = .75 * (1 - z.^2); % epanechnikov-like weights
k = k ./ sum(k);
Z = filter2(k'*k,Y);

function varargout = gsp(x,y,c,ms,ColorMap)
%Graphs scattered poits
map = ColorMap;
ind = fix((c-min(c))/(max(max(c,2))-min(c))*(size(map,1)-1))+1;
% ind = fix((c-min(c))/(max(c)-min(c))*(size(map,1)-1))+1;
h = [];
%much more efficient than matlab's scatter plot
for k=1:size(map,1)
    if any(ind==k)
%         h(end+1) = line('Xdata',x(ind==k),'Ydata',y(ind==k), ...
%             'LineStyle','none','Color',map(k,:), ...
%             'Marker','.','MarkerSize',ms);
        h(end+1) = line('Xdata',x(ind==k),'Ydata',y(ind==k), ...
            'LineStyle','none','Color',map(k,:), ...
            'Marker','o','MarkerSize',ms,'MarkerFaceColor',map(k,:));
    end
end
if nargout==1
    varargout{1} = h;
end
return

