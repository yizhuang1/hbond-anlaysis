#!/usr/bin/env python
import sys,os,fnmatch,itertools,pickle
import os.path,datetime,time
from glob import glob
import numpy as np
from random import *

#removed matplotlib here for keeneland, but it is not
#needed in general...
#import matplotlib
#import matplotlib.pyplot as plt

my_dir = os.path.abspath(os.path.dirname(__file__))
num=sys.argv[1]

def pack_pkl(stage):
    ''' Combine bonding pickles into 1 file.
        hb_protein-protein.pkl  =>  01-sd_hb.pkl
    '''
    with open('%s-sd_hb.pkl' % stage,'wb') as hpk:
        for path in glob(os.path.join(my_dir,'%s/*/*-hb_pr*pr*.pkl.*' % stage)):
        #path to each hb.pkl in the sub directories
        #print 'path',path
            dct_s_hb={}
            #print path
            seed = path.split('.')[-1]
            sample_i = pickle.load(open(path,'rb'))
            if len(sample_i)==100:
            #print 'sample_i',sample_i
            #QUESTION: what changes the variable there; didn't realize it was a variable
                dct_s_hb[seed]=[sample_i]
                #os.remove(path) #HAILEY 6-12-14
                #print 'dct_s_hb',dct_s_hb[seed]
                #print 'hpk',hpk
            if len(dct_s_hb)>0:
                pickle.dump(dct_s_hb,hpk,protocol=2) 

    ''' Combine bonding pickles into 1 file.
        protein-water           =>  01-sd_wp.pkl
    '''
    if my_dir.split('/')[-2].split('.')[1]=='exp':
        with open('%s-sd_wp.pkl' % stage,'wb') as wpk:
            for path in glob(os.path.join(my_dir,'%s/*/*-hb_pr*resname*.pkl.*' % stage)):
                dct_s_whb={}
                print(path)
                seed = path.split('.')[-1]
                sample_i = pickle.load(open(path,'rb'))
                if len(sample_i)==100:
                    dct_s_whb[seed]=[sample_i]
                    #os.remove(path) #HAILEY 6-12-14
                if len(dct_s_hb)>0:
                    pickle.dump(dct_s_whb,wpk,protocol=2)

# main call
pack_pkl(num)
