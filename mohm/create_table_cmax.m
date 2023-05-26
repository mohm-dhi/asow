
clear
clc


for i =1:7
    path = '\\\\usden1-stor\\Projects\\41806529\\08_Results\\_met_results\\03 Extreme\\3.4 Waves\\3.4.3 Hmax and Cmax\\ASOW%s\\1000\\Cmax\\Atlantic_Shores_Offshore_Wind_Forristall_Cmax_%s.mat';
msl_dir = sprintf(path, num2str(i), 'MSL_Directional');
swl_dir = sprintf(path, num2str(i), 'SWL_Directional');
msl_month = sprintf(path, num2str(i), 'MSL_Monthly');
swl_month = sprintf(path, num2str(i), 'SWL_Monthly');

msl_dir_data = load(msl_dir);
name = fieldnames(msl_dir_data);
EVA = getfield(msl_dir_data, name{1});
filename = replace(msl_dir, '.mat', '.xlsx');
[filepath,name,ext] = fileparts(filename);
for t = 1:length(EVA.T), ylabs(t) = {num2str(EVA.T(t))}; end
subtitle{1} = 'All';
cd(filepath);
m_table(['MSL_Directional_' EVA.filename],{['P' num2str(i)] subtitle{1}},{'T_R [years]'},{[EVA.label(1) '_{max}' '  (' EVA.unit ')']},ylabs,EVA.title([max(1,EVA.rep_dirs(end)-length(EVA.rep_dirs)) EVA.rep_dirs(1:end-1)+1]),circshift(EVA.R_mx(:,1:end,2),[0 1]),EVA.ascii,'Precision',EVA.ndec); 


swl_dir_data = load(swl_dir);
name = fieldnames(swl_dir_data);
EVA = getfield(swl_dir_data, name{1});
filename = replace(swl_dir, '.mat', '.xlsx');
[filepath,name,ext] = fileparts(filename);
for t = 1:length(EVA.T), ylabs(t) = {num2str(EVA.T(t))}; end
subtitle{1} = 'All';
cd(filepath);
m_table(['SWL_Directional_' EVA.filename],{['P' num2str(i)] subtitle{1}},{'T_R [years]'},{[EVA.label(1) '_{max}' '  (' EVA.unit ')']},ylabs,EVA.title([max(1,EVA.rep_dirs(end)-length(EVA.rep_dirs)) EVA.rep_dirs(1:end-1)+1]),circshift(EVA.R_mx(:,1:end,2),[0 1]),EVA.ascii,'Precision',EVA.ndec); 



msl_month_data = load(msl_month);
name = fieldnames(msl_month_data);
EVA = getfield(msl_month_data, name{1});
filename = replace(msl_month, '.mat', '.xlsx');
[filepath,name,ext] = fileparts(filename);
for t = 1:length(EVA.T), ylabs(t) = {num2str(EVA.T(t))}; end
subtitle{1} = 'All';
cd(filepath);
m_table(['MSL_Monthly_' EVA.filename],{['P' num2str(i)] subtitle{1}},{'T_R [years]'},{[EVA.label(1) '_{max}' '  (' EVA.unit ')']},ylabs,EVA.title([max(1,EVA.rep_dirs(end)-length(EVA.rep_dirs)) EVA.rep_dirs(1:end-1)+1]),circshift(EVA.R_mx(:,1:end,2),[0 1]),EVA.ascii,'Precision',EVA.ndec); 


swl_month_data = load(swl_month);
name = fieldnames(swl_month_data);
EVA = getfield(swl_month_data, name{1});
filename = replace(swl_month, '.mat', '.xlsx');
[filepath,name,ext] = fileparts(filename);
for t = 1:length(EVA.T), ylabs(t) = {num2str(EVA.T(t))}; end
subtitle{1} = 'All';
cd(filepath);
m_table(['SWL_Monthly_' EVA.filename],{['P' num2str(i)] subtitle{1}},{'T_R [years]'},{[EVA.label(1) '_{max}' '  (' EVA.unit ')']},ylabs,EVA.title([max(1,EVA.rep_dirs(end)-length(EVA.rep_dirs)) EVA.rep_dirs(1:end-1)+1]),circshift(EVA.R_mx(:,1:end,2),[0 1]),EVA.ascii,'Precision',EVA.ndec); 




end