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

dst = r'z:\03 Operational Conditions\3.4 Waves\3.4.1 Annual and Monthly Statistics\P%d'

for stat in ['Directional','Monthly']:

    for k in range(1,8):

        src = r'\\USDEN1-STOR.DHI.DK\Projects\41806529\08_Results\_met_results\02 Operational\2.4 Waves\2.4.1 Annual and Monthly Directional\%sStats\P%d' % (stat, k)
        
        # %% plots
        for part in ['Hm0_Total','Hm0_Sea','Hm0_Swell']:
            
        
            src_k = os.path.join(src, part)
            dst_k = os.path.join(dst % k, part)
            
            if not os.path.isdir(dst_k):
                os.makedirs(dst_k)
                
            # list and copy files
            for src_file in glob.glob(os.path.join(src_k,'*.jpg')):
                print(src_file)
                dst_file = os.path.basename(src_file.replace('ASOW','P')).replace('.','').replace('jpg','.jpg').replace('_rose','').replace('Scatter_','Rose_')
                dst_file = os.path.join(dst_k, dst_file)
                shutil.copyfile(src_file,dst_file)
              
    
        # %% tables
        for part in ['Hm0_Total','Hm0_Sea','Hm0_Swell']:
    
            src_k = os.path.join(src, part)
            dst_k = os.path.join(dst % k, part)
            
            if not os.path.isdir(dst_k):
                os.mkdir(dst_k)
                
            # list and copy files
            for src_file in glob.glob(os.path.join(src_k,'*.xls')):
                print(src_file)
                dst_file = os.path.basename(src_file.replace('ASOW','P')).replace('.xls','_%s.xls' % stat)
                dst_file = os.path.join(dst_k, dst_file)
                shutil.copyfile(src_file,dst_file)        
            