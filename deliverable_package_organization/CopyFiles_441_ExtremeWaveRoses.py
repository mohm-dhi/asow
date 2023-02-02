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

dst = r'z:\04 Extreme\4.4 Waves\4.4.1 Peak Wave Roses\P%d'

for k in range(1,8):
    
    src = r'\\USDEN1-STOR.DHI.DK\Projects\41806529\08_Results\_met_results\03 Extreme\3.4 Waves\3.4.1 Peak Wave Rose\P%d' % (k)
    
    src_k = os.path.join(src)
    dst_k = os.path.join(dst % k)
    
    if not os.path.isdir(dst_k):
        os.makedirs(dst_k)
        

    # list and copy files - figures
    for src_file in glob.glob(os.path.join(src_k,'*Omni_Rose.jpg')):
        print(src_file)
        dst_file = os.path.basename(src_file.replace('ASOW','P')).replace('.','').replace('jpg','.jpg').replace('_Data','').replace(')directional',')')
        dst_file = os.path.join(dst_k, dst_file)
        shutil.copyfile(src_file,dst_file)
        
