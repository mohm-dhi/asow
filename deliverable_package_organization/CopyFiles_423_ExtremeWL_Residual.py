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

dst = r'z:\04 Extreme\4.2 Water Levels\4.2.3 Residual\P%d\%s'

for k in range(1,8):
    
    for part in ['HWL','LWL']:
    
        src = r'\\USDEN1-STOR.DHI.DK\Projects\41806529\08_Results\_met_results\03 Extreme\3.2 Water Levels\3.2.2 Residual\P%d\%s' % (k, part)
        
        src_k = os.path.join(src)
        dst_k = os.path.join(dst % (k,part))
        
        if not os.path.isdir(dst_k):
            os.makedirs(dst_k)
            
        # list and copy files - figures
        for src_file in glob.glob(os.path.join(src_k,'*.jpg')):
            print(src_file)
            dst_file = os.path.basename(src_file.replace('ASOW','P')).replace('.','').replace('jpg','.jpg')
            dst_file = os.path.join(dst_k, dst_file)
            shutil.copyfile(src_file,dst_file)
            
        # list and copy files - tables
        for src_file in glob.glob(os.path.join(src_k,'*.txt')):
            print(src_file)
            dst_file = os.path.basename(src_file.replace('ASOW','P'))
            dst_file = os.path.join(dst_k, dst_file)
            shutil.copyfile(src_file,dst_file)
