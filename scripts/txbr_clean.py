#!/usr/bin/python

import os

from txbr import ALIGN_DIR, \
                 FILTER_DIR, \
                 SETUP_DIR, \
                 BACKPROJECTION_DIR, \
                 CONTOUR_DIR

def rm_rf(d):
    '''A routine that allows to remove everything (directory and files) recursively'''
    
    for path in (os.path.join(d,f) for f in os.listdir(d)):
        if os.path.isdir(path):
            rm_rf(path)
        else:
            os.unlink(path)
            
    os.rmdir(d)


def main(directory):
    
    var = raw_input("Are you sure? (Y or N): ")
    if var!='Y': 
        return
    
    files = os.listdir(directory)
    
    for file in files:
        
        if file.endswith('.st'): continue
        if file.endswith('.preali'): continue
        if file.endswith('.fid'): continue
        if file.endswith('.rawtlt'): continue
        if file.endswith('.prexg'): continue
        if file.endswith('.mosaic'): continue
        if file.endswith('.common'): continue
        if file.endswith('.xg'): continue
        if file.endswith('.prexf'): continue
        if file.endswith('.xf'): continue
        if file.endswith('.contmod'): continue
        if file.endswith('.pbs'): continue
        if file.endswith('.py'): continue
        if file.endswith('.cfg'): continue
        
        if os.path.isdir(file)==True: 
            if file==ALIGN_DIR or file==FILTER_DIR or file==SETUP_DIR or file==BACKPROJECTION_DIR: 
                rm_rf(file)
            continue
        
        os.remove(file)
        
main('.')