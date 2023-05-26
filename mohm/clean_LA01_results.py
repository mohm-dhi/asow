import os
from os import listdir
from os.path import join

path = r"\\usden1-nas3\41806593\07-timeseries\JEVA\03-Output\41806593\ASOW_FS3\03_ReturnValues\diagnostics_plots"

LA02_files = [join(path, f) for f in listdir(path) if 'LA_02' in f]
LA01_files = [join(path, f) for f in listdir(path) if 'LA_01' in f]

diff = [f for f in LA01_files if f.replace('LA_01', 'LA_02') not in LA02_files]

print(len(diff))
for file in diff:
    os.remove(file)