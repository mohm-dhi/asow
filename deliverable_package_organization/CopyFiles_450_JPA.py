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

dst = r'z:\04 Extreme\4.5 Joint Probability\4.5.%d %s\P%d'

for n, part1 in enumerate(['Wind SPD vs Hm0','Hm0 vs Water Level','Hm0 vs Depth Avg Current Speed']):
    
    n += 1
    
    if 'Current' in part1:
        part2 = 'Hm0 vs Currents'
    else:
        part2 = part1
    

    for k in range(1,8):
        
        src = r'\\USDEN1-STOR.DHI.DK\Projects\41806529\08_Results\_met_results\03 Extreme\3.5 JPA\3.5.%d %s\P%d' % (n,part2,k)
        
        src_k = os.path.join(src)
        dst_k = os.path.join(dst % (n,part1,k))
        
        if not os.path.isdir(dst_k):
            os.makedirs(dst_k)            
    
        # list and copy files - figures
        for src_file in glob.glob(os.path.join(src_k,'*_from_plot1.png')):
            print(src_file)
            dst_file = os.path.basename(src_file.replace('ASOW','P')).replace('.','').replace('png','.png').replace('_from_plot1','')
            dst_file = os.path.join(dst_k, dst_file)
            shutil.copyfile(src_file,dst_file)
        
        # # list and copy files - tables
        # for src_file in glob.glob(os.path.join(src_k,'*.csv')):
            
        #     print(src_file)
        #     dst_file = os.path.basename(src_file.replace('ASOW','P'))
        #     dst_file = os.path.join(dst_k, dst_file)
        #     shutil.copyfile(src_file,dst_file)
        
