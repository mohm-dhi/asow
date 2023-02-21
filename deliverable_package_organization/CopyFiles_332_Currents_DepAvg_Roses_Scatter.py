# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 13:11:01 2023

Copy files

- personalized for each type of plots

@author: alka
"""

import os
import shutil
import glob

# %% Wind Normal Stats - Hub Height

src = r'\\USDEN1-STOR.DHI.DK\Projects\41806529\08_Results\_met_results\02 Operational\2.3 Currents\2.3.1 Annual and Monthly Depth Avg\P%d'
dst = r'z:\03 Operational Conditions\3.3 Currents\3.3.2 Depth Avg Currents - Roses and Frequency Tables\P%d'

for k in range(1,8):
    
    # %% roses
    for part in ['Residual','Tide','Total']:
    
        src_k = os.path.join(src % k,'Rose',part)
        dst_k = os.path.join(dst % k, part)
        
        if not os.path.isdir(dst_k):
            os.makedirs(dst_k)
            
        # list and copy files
        for src_file in glob.glob(os.path.join(src_k,'*Depth_Averaged*.jpg')):
            print(src_file)
            dst_file = os.path.basename(src_file.replace('ASOW','P')).replace('.','').replace('jpg','.jpg').replace('_rose','').replace('Scatter_','Rose_')
            dst_file = os.path.join(dst_k, dst_file)
            shutil.copyfile(src_file,dst_file)
          

    # %% scatter tables
    for part in ['Residual','Tide','Total']:

        src_k = os.path.join(src % k,'Table',part)
        dst_k = os.path.join(dst % k, part)
        
        if not os.path.isdir(dst_k):
            os.mkdir(dst_k)
            
        # list and copy files
        for src_file in glob.glob(os.path.join(src_k,'*Scatter*Depth_Averaged*.txt')):
            print(src_file)
            dst_file = os.path.basename(src_file.replace('ASOW','P')).replace('.','').replace('txt','.txt').replace('_rose','')
            dst_file = os.path.join(dst_k, dst_file)
            shutil.copyfile(src_file,dst_file)        
            