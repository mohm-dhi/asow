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

dst = r'z:\04 Extreme\4.3 Currents\4.3.2 Total Current\P%d'


depin = ['0. Depth Averaged',
        '1. Near Surface',
        '2. 75',
        '3. Mid',
        '4. 25',
        '5. Near Seabed']

depou = ['DepthAveraged',
        'Surface',
        '75perc',
        'Mid_Column',
        '25perc',
        'Near_Bed']

for d1,d2 in zip(depin,depou):
    
    for part in ['Directional','Monthly Omnidirectional']:
    
        for k in range(1,8):
            
            src = r'\\USDEN1-STOR.DHI.DK\Projects\41806529\08_Results\_met_results\03 Extreme\3.3 Currents\3.3.2 Total\P%d\%s\%s' % (k,d1,part)
            
            src_k = os.path.join(src)
            dst_k = os.path.join(dst % k, d2)
            
            if not os.path.isdir(dst_k):
                os.makedirs(dst_k)            
        
            # list and copy files - figures
            for src_file in glob.glob(os.path.join(src_k,'*.jpg')):
                print(src_file)
                dst_file = os.path.basename(src_file.replace('ASOW','P')).replace('.','').replace('jpg','.jpg').replace('Directionaldirectional','Directional').replace('Monthlydirectional','Monthly')
                dst_file = os.path.join(dst_k, dst_file)
                shutil.copyfile(src_file,dst_file)
                
            # list and copy files - tables
            for src_file in glob.glob(os.path.join(src_k,'*.txt')):
                
                if 'Monthlyomni' in src_file:
                    continue
                
                print(src_file)
                dst_file = os.path.basename(src_file.replace('ASOW','P').replace('Directionaldirectional','Directional').replace('Monthlydirectional','Monthly'))
                dst_file = os.path.join(dst_k, dst_file)
                shutil.copyfile(src_file,dst_file)
        
