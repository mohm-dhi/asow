# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 09:29:33 2023

@author: alka
"""

import os

startpath = r'z:/'

def list_files(startpath, fname):
    
    with open(fname,'w') as fid:
    
        for root, dirs, files in os.walk(startpath):
            
            if 'Scripts' in root:
                continue
            
            level = root.replace(startpath, '').count(os.sep)
            indent = ' ' * 4 * (level)
            fid.write('{}{}/'.format(indent, os.path.basename(root)))
            fid.write('\n')
            subindent = ' ' * 4 * (level + 1)
            for f in files:
                
                if 'Thumbs.db' in f:
                    continue

                fid.write('{}{}'.format(subindent, f))
                fid.write('\n')
            
list_files(startpath,'FolderContent.txt')