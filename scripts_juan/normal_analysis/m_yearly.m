function [X] = m_yearly(X,year_arr)

% Copyright DHI ©, Denmark, All rights reserved
% 2010-11-09 - Patrick Dich Grode, pdg@dhigroup.com

% Add titles, colors and line types
% X.title  = [X.title  {'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec'}];
% X.color  = [X.color  ;X.ColorOrder(2,:);   X.ColorOrder(3,:);   X.ColorOrder(4,:);   X.ColorOrder(5,:);   X.ColorOrder(2,:);   X.ColorOrder(3,:);   X.ColorOrder(4,:);   X.ColorOrder(5,:);   X.ColorOrder(2,:);   X.ColorOrder(3,:);   X.ColorOrder(4,:);   X.ColorOrder(5,:);];
% X.line   = [X.line   {'-'   '-'   '-'   '-'   '--'  '--'  '--'  '--'  ':'   ':'   ':'   ':'}];
% X.marker = [X.marker {'o'   's'   '^'   'x'   's'   '^'   'x'   'o'   '^'   'x'   'o'   's' }];

% Add tag and data for each new column
%  keyboard
for i = 1:length(year_arr)
    % Add group
    X.group = [X.group 'Yearly'];
    % Initialize column with nan values
    X.data(:,end+1) = nan;
    % Find times within limits
    t_vec = datevec(X.time);
    % Find current year
    r = find(t_vec(:,1) == year_arr(i));
    % Change the values of the new column
    X.data(r,end) = X.data(r);
    % Add year yo title cell
    X.title  = [X.title  {num2str(year_arr(i))}];
end 