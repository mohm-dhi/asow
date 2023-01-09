function [sw_arr] = read_mood_sw_csv(fdir,fname,headerlines)

% fdir = 'C:\Users\jngz\OneDrive - DHI\2022\41806529-AtlanticShores\Data\Points to be delivered SW_HD\';
% fname = 'P1_US_EastCoast_HD_-73.945_39.651_24.6_3051.4_1979-01-01_2021-12-31_.csv';
% headerlines = 16;

% read MOOD csv into a table, time in datetime
tt = readtable([fdir fname],'NumHeaderLines',headerlines);

% convert datetime to datenum, split wl and currents/dir
sw_arr = [datenum(tt{:,1}) tt{:,2:end}];

end