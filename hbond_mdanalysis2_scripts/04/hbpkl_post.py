#!/usr/bin/env python
from sys import argv
import numpy as np
import os,sys,pickle
import os.path
from glob import glob


my_dir = os.path.abspath(os.path.dirname(__file__))
for path in glob(os.path.join(my_dir,'daOut.dcd.*')):
    file=path.split('/')[-1]
    withoutseed=file[0:9]
    seed=file[10:15]
    os.rename(file,withoutseed)
    for path2 in glob(os.path.join(my_dir,'*-tef.dat.%s' %(seed))):
        number = int(path2.split('/')[-1].split('-')[0])
        os.system('python3 ../hb.py %d %s' % (number, seed))
    os.rename(withoutseed,file)
