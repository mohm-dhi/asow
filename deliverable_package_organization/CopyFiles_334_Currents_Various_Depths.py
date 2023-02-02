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

for depth in ['Surface','Near_Bed','Mid_Column']:

    for k in range(1,8):
        
        src = r'\\USDEN1-STOR.DHI.DK\Projects\41806529\08_Results\_met_results\02 Operational\2.3 Currents\2.3.1 Annual and Monthly Depth Avg\P%d'  % k
        dst = r'z:\03 Operational Conditions\3.3 Currents\3.3.4 Currents at Various Depths - Statistics\P%d\%s' % (k, depth)
              
        # %% probability
        for part in ['Residual','Tide','Total']:
    
            src_k = os.path.join(src,'Probability', part)
            dst_k = os.path.join(dst, part, 'Monthly_Probability')
            
            if not os.path.isdir(dst_k):
                os.makedirs(dst_k)
                
            # list and copy files
            for src_file in glob.glob(os.path.join(src_k,'*%s*.jpg' % depth.replace('_',' '))):
                print(src_file)
                dst_file = os.path.basename(src_file.replace('ASOW','P')).replace('.','').replace('jpg','.jpg')
                dst_file = os.path.join(dst_k, dst_file)
                shutil.copyfile(src_file,dst_file)
                     
                
        # %% general stats
        for part,partname in zip(['Residual','Tide','Total'],
                     ['Res','Tide','Total']):
    
            src_k = os.path.join(src)
            dst_k = os.path.join(dst, part)
            
            if not os.path.isdir(dst_k):
                os.mkdir(dst_k)
                
            # list and copy files
            for src_file in glob.glob(os.path.join(src_k,'*%s*%s*STATS.xlsx' % (depth.replace('_',' '),partname))):
                print(src_file)
                dst_file = os.path.basename(src_file.replace('ASOW%d' % k,'P%d_' % k)).replace('.','').replace('xlsx','.xlsx')
                dst_file = os.path.join(dst_k, dst_file)
                shutil.copyfile(src_file,dst_file)       

            # %% roses
            src_k = os.path.join(src,'Rose',part)
            # list and copy files
            for src_file in glob.glob(os.path.join(src_k,'*%s*.jpg' % depth)):
                print(src_file)
                dst_file = os.path.basename(src_file.replace('ASOW','P')).replace('.','').replace('jpg','.jpg').replace('_rose','').replace('Scatter_','Rose_')
                dst_file = os.path.join(dst_k, dst_file)
                shutil.copyfile(src_file,dst_file)
              
    
            # %% scatter tables
            src_k = os.path.join(src,'Table',part)
            # list and copy files
            for src_file in glob.glob(os.path.join(src_k,'*Scatter*%s*.txt' % depth)):
                print(src_file)
                dst_file = os.path.basename(src_file.replace('ASOW','P')).replace('.','').replace('txt','.txt').replace('_rose','')
                dst_file = os.path.join(dst_k, dst_file)
                shutil.copyfile(src_file,dst_file)                    