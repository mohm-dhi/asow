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
import pandas as pd

dst = r'z:\04 Extreme\4.5 Joint Probability\4.5.%d %s\P%d'

for n, part1 in enumerate(['Wind SPD vs Hm0','Hm0 vs Water Level','Hm0 vs Depth Avg Current Speed']):
    
    n += 1
    
    if 'Current' in part1:
        part2 = 'Hm0 vs Currents'
    else:
        part2 = part1
        
    # define label
    if part1 == 'Wind SPD vs Hm0':
        lab1 = 'WS_{HubHeight} [m/s]'
        lab2 = 'H_{m0} [m]'
    elif part1 == 'Hm0 vs Water Level':
        lab1 = 'H_{m0} [m]'
        lab2 = 'WL_{Total} [m]'
    elif part1 == 'Hm0 vs Depth Avg Current Speed':
        lab1 = 'H_{m0} [m]'
        lab2 = 'CS_{Total} [m/s]'    
    else:
        aaaa

    for k in range(1,8):
        
        src = r'\\USDEN1-STOR.DHI.DK\Projects\41806529\08_Results\_met_results\03 Extreme\3.5 JPA\3.5.%d %s\P%d' % (n,part2,k)
        
        src_k = os.path.join(src)
        dst_k = os.path.join(dst % (n,part1,k))
        
        if not os.path.isdir(dst_k):
            os.makedirs(dst_k)            
        
        # list and copy files - tables
        for src_file in glob.glob(os.path.join(src_k,'*.csv')):
            
            df = pd.read_csv(src_file)
            df.columns = ['T_R [years]',
                          lab1.split()[0],
                          '%s 5%%' % lab2.split()[0],
                          '%s 50%%' % lab2.split()[0],
                          '%s 95%%' % lab2.split()[0]];
            
            print(src_file)
            dst_file = 'P%d_' % k + os.path.basename(src_file)
            dst_file = os.path.join(dst_k, dst_file.replace('.csv','.txt'))
            
            with open(dst_file,'w') as f:
                f.write('P%d\n' % k)
                f.write('Joint Probability Analysis - %s x %s\n\n' % (lab1,lab2))
                
                f.write('%11s\t%14s\t%14s\t%14s\t%14s\n' % tuple(df.columns))
                
                for i, row in df.iterrows():
                    f.write('%11.0f\t%14.1f\t%14.2f\t%14.2f\t%14.2f\n' % tuple(row.values))

            
            
        
