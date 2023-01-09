clear
clc

file = "C:\Users\mohm\OneDrive - DHI\Documents\Project Files\41806529 - Atlantic Shores\Mode_Obs_Comparison\Data\Processed_Data\Loc4\QC_datalistings_aqd_SV3_WS201.csv"
out_file = strrep(file,'.csv', '_depth_averaged.csv');

data = readtable(file);
head(data)
values = table2array(data(:, 2:end));
uv = values*0; 

for i = 1:size(values, 2)/2
    [u, v] = spddir2uv(values(:, 2*i -1), values(:, 2* i ) );
    uv(:, i) = u;
    uv(:, i + size(values, 2)/2) = v;
    
end


u_avg = nanmean(uv(:, 1:size(values, 2)/2), 2);
v_avg = nanmean(uv(:, size(values,2)/2 + 1 : end) , 2);

[spd, dir] = uv2spddir(u_avg, v_avg, 'to');


d = data(:, 1);
d.Uavg = u_avg;
d.Vavg = v_avg;
d.spd = spd;
d.dir = dir;

writetable(d, out_file)

