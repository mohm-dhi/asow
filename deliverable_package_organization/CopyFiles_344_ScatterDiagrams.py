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

dst = r'z:\03 Operational Conditions\3.4 Waves\3.4.4 Scatter Diagrams\%s\P%d'

# %% Hm0_T02

for k in range(1,8):
    
    # Hm0_T02
    src = r'\\USDEN1-STOR.DHI.DK\Projects\41806529\08_Results\_met_results\02 Operational\2.4 Waves\2.4.4 Scatter Diagrams\P%d\Hm0_T02\Directional\Hm0_Total' % k
    
    src_k = os.path.join(src)
    dst_k = os.path.join(dst % ('Hm0_T02', k))
    
    if not os.path.isdir(dst_k):
        os.makedirs(dst_k)
        
    # list and copy files - figures
    for src_file in glob.glob(os.path.join(src_k,'*.jpg')):
        print(src_file)
        dst_file = os.path.basename(src_file.replace('ASOW','P')).replace('.','').replace('jpg','.jpg')
        dst_file = os.path.join(dst_k, dst_file)
        shutil.copyfile(src_file,dst_file)
        
    # list and copy files - tables
    for src_file in glob.glob(os.path.join(src_k,'*.xlsx')):
        print(src_file)
        dst_file = os.path.basename(src_file.replace('ASOW','P'))
        dst_file = os.path.join(dst_k, dst_file)
        shutil.copyfile(src_file,dst_file)        
        
# %% Hm0_WSHub

for k in range(1,8):
    
    src = r'\\USDEN1-STOR.DHI.DK\Projects\41806529\08_Results\_met_results\02 Operational\2.4 Waves\2.4.4 Scatter Diagrams\P%d\Hm0_WSHub\Directional' % k
    
    for part in ['Hm0_Total','Hm0_Sea','Hm0_Swell']:
        src_k = os.path.join(src, part)
        dst_k = os.path.join(dst % ('Hm0_WSHub',k), part)
        
        if not os.path.isdir(dst_k):
            os.makedirs(dst_k)
            
        # list and copy files - figures
        for src_file in glob.glob(os.path.join(src_k,'*.jpg')):
            print(src_file)
            dst_file = os.path.basename(src_file.replace('ASOW','P')).replace('.','').replace('jpg','.jpg')
            dst_file = os.path.join(dst_k, dst_file)
            shutil.copyfile(src_file,dst_file)
            
        # list and copy files - tables
        for src_file in glob.glob(os.path.join(src_k,'*.xlsx')):
            print(src_file)
            dst_file = os.path.basename(src_file.replace('ASOW','P'))
            dst_file = os.path.join(dst_k, dst_file)
            shutil.copyfile(src_file,dst_file)                
            
# %% Hm0_Tp

for k in range(1,8):
    
    src = r'\\USDEN1-STOR.DHI.DK\Projects\41806529\08_Results\_met_results\02 Operational\2.4 Waves\2.4.4 Scatter Diagrams\P%d\Hm0_Tp\Directional' % k
    
    for part in ['Hm0_Total','Hm0_Sea','Hm0_Swell']:
        src_k = os.path.join(src, part)
        dst_k = os.path.join(dst % ('Hm0_Tp',k), part)
        
        if not os.path.isdir(dst_k):
            os.makedirs(dst_k)
            
        # list and copy files - figures
        for src_file in glob.glob(os.path.join(src_k,'*.jpg')):
            print(src_file)
            dst_file = os.path.basename(src_file.replace('ASOW','P')).replace('.','').replace('jpg','.jpg')
            dst_file = os.path.join(dst_k, dst_file)
            shutil.copyfile(src_file,dst_file)
            
        # list and copy files - tables
        for src_file in glob.glob(os.path.join(src_k,'*.xlsx')):
            print(src_file)
            dst_file = os.path.basename(src_file.replace('ASOW','P'))
            dst_file = os.path.join(dst_k, dst_file)
            shutil.copyfile(src_file,dst_file)                
            
# %% WDir_MWD

for k in range(1,8):
    
    src = r'\\USDEN1-STOR.DHI.DK\Projects\41806529\08_Results\_met_results\02 Operational\2.4 Waves\2.4.4 Scatter Diagrams\P%d\WDir_MWD\Monthly' % k
    
    for part in ['MWD_Total','MWD_Sea','MWD_Swell']:
        src_k = os.path.join(src, part)
        dst_k = os.path.join(dst % ('WDir_MWD', k), part)
        
        if not os.path.isdir(dst_k):
            os.makedirs(dst_k)
            
        # list and copy files - figures
        for src_file in glob.glob(os.path.join(src_k,'*.jpg')):
            print(src_file)
            dst_file = os.path.basename(src_file.replace('ASOW','P')).replace('.','').replace('jpg','.jpg')
            dst_file = os.path.join(dst_k, dst_file)
            shutil.copyfile(src_file,dst_file)
            
        # list and copy files - tables
        for src_file in glob.glob(os.path.join(src_k,'*.xlsx')):
            print(src_file)
            dst_file = os.path.basename(src_file.replace('ASOW','P'))
            dst_file = os.path.join(dst_k, dst_file)
            shutil.copyfile(src_file,dst_file)                            