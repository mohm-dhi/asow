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

dst = r'z:\06 Env Parameters\6.1 Seawater and Meteorological Parameters\6.1.%d %s'

parts = ['Air Density',
        'Air Temp',
        'Humidity',
        'Salinity',
        'Seawater Density',
        'Seawater Temp',
        'Solar Rad']

order  = [2,
          1,
          3,
          7,
          5,
          6,
          4]


for n, part in zip(order,parts):

    for k in range(1,8):
        
        src = r'\\USDEN1-STOR.DHI.DK\Projects\41806529\08_Results\_met_results\05 Env Parameters\5.1 Met Stats\%s' % (part)
        
        src_k = os.path.join(src)
        dst_k = os.path.join(dst % (n,part))
        
        if not os.path.isdir(dst_k):
            os.makedirs(dst_k)            
    
        # list and copy files - figures
        for src_file in glob.glob(os.path.join(src_k,'*.jpg')):
            print(src_file)
            dst_file = os.path.basename(src_file.replace('ASOW3','CFSR')).replace('.','').replace('jpg','.jpg')
            dst_file = os.path.join(dst_k, dst_file)
            shutil.copyfile(src_file,dst_file)

        # list and copy files - figures
        for src_file in glob.glob(os.path.join(src_k,'*.txt')):
            print(src_file)
            dst_file = os.path.basename(src_file.replace('ASOW3','CFSR')).replace('.','').replace('txt','.txt')
            dst_file = os.path.join(dst_k, dst_file)
            shutil.copyfile(src_file,dst_file)
        
        # # list and copy files - tables
        # for src_file in glob.glob(os.path.join(src_k,'*.csv')):
            
        #     print(src_file)
        #     dst_file = os.path.basename(src_file.replace('ASOW','P'))
        #     dst_file = os.path.join(dst_k, dst_file)
        #     shutil.copyfile(src_file,dst_file)
        
