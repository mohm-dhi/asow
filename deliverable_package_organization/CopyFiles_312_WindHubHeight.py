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

src = r'\\USDEN1-STOR.DHI.DK\Projects\41806529\08_Results\_met_results\02 Operational\2.1 Wind\2.1.1 Annual and Monthly\P%d\HubHeight\FullDataSet'
dst = r'\\USDEN1-STOR.DHI.DK\Projects\41806529\08_Results\Client_Deliverable\v001\03 Operational Conditions\3.1 Wind\3.1.2 Stats - Hub Height\P%d'

for k in range(1,8):
    
    # wind roses
    src_k = os.path.join(src % k,'Rose')
    dst_k = os.path.join(dst % k,'Roses')
    
    if not os.path.isdir(dst_k):
        os.mkdir(dst_k)
        
    # list and copy files
    for src_file in glob.glob(os.path.join(src_k,'*.jpg')):
        print(src_file)
        dst_file = os.path.basename(src_file.replace('ASOW','P'))
        dst_file = os.path.join(dst_k, dst_file)
        shutil.copyfile(src_file,dst_file)
        
    # directional probability
    src_k = os.path.join(src % k,'Probabilty')
    dst_k = os.path.join(dst % k,'Directional_Probability')
    
    if not os.path.isdir(dst_k):
        os.mkdir(dst_k)
        
    # list and copy files
    for src_file in glob.glob(os.path.join(src_k,'*.jpg')): 
        if src_file[-5].isdigit():
            dst_file = os.path.basename(src_file.replace('ASOW','P'))
            dst_file = os.path.join(dst_k, dst_file)
            shutil.copyfile(src_file,dst_file)
    
    # monthly probability
    src_k = os.path.join(src % k,'Probabilty')
    dst_k = os.path.join(dst % k,'Monthly_Probability')
    
    if not os.path.isdir(dst_k):
        os.mkdir(dst_k)
        
    # list and copy files
    for src_file in glob.glob(os.path.join(src_k,'*.jpg')):
        if not src_file[-5].isdigit():
            dst_file = os.path.basename(src_file.replace('ASOW','P'))
            dst_file = os.path.join(dst_k, dst_file)
            shutil.copyfile(src_file,dst_file)
    
    # scatter tables
    src_k = os.path.join(src % k,'Table')
    dst_k = os.path.join(dst % k,'Scatter_Tables')
    
    if not os.path.isdir(dst_k):
        os.mkdir(dst_k)
        
    # list and copy files
    for src_file in glob.glob(os.path.join(src_k,'*Scatter*.xlsx')):
        if not src_file[-5].isdigit():
            dst_file = os.path.basename(src_file.replace('ASOW','P'))
            dst_file = os.path.join(dst_k, dst_file)
            shutil.copyfile(src_file,dst_file)    
    
    # general stats
    src_k = os.path.join(src % k,'..')
    dst_k = os.path.join(dst % k,'General_Statistics')
    
    if not os.path.isdir(dst_k):
        os.mkdir(dst_k)
        
    # list and copy files
    for src_file in glob.glob(os.path.join(src_k,'*STATS*.xlsx')):
        if not src_file[-5].isdigit():
            dst_file = os.path.basename(src_file.replace('ASOW','P'))
            dst_file = os.path.join(dst_k, dst_file)
            shutil.copyfile(src_file,dst_file)    
