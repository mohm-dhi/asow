file = "\\USDEN1-STOR.DHI.DK\Projects\41806529\01_Client_supplied_data\3 Data\Data\Loc6a\144041-QC_datalistings_SV1_03\144041-QC_datalistings_spectra_SV1_03.LIS"

out_file = strrep(file,'.LIS', '.csv')
data = readtable(file, 'FileType','text', 'Range', 'B23');
head(data)
data.Properties.VariableNames{1} = 'Time';
data.Time = datetime(string(data{:,1}), 'format', 'yyyy-MM-dd HH:mm:ss');
writetable(data, out_file)

head(data)
