#!/usr/bin/env python
import sys,os,fnmatch,itertools,pickle,time
import os.path
from glob import glob
import numpy as np
from random import *
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import LSQUnivariateSpline
import plottools
from matplotlib.ticker import MultipleLocator

my_dir = os.path.abspath(os.path.dirname(__file__))
fig,ax = plt.subplots()
plt.tick_params(axis='both',bottom=True,top=True,left=True,right=True,direction='in',which='both',grid_color='blue')
class mdict(dict):
    def __setitem__(self,key,value):
        self.setdefault(key,[]).append(value)
def print_dict(dct):
    for key,val in dct.items():
        print key,val
        print ''
    return key
count = 0
for path in glob(os.path.join(my_dir,'*-sfwf.pkl*')):
    count +=1
num = str(count).zfill(2)
if '00' == num:
    print num,'no data acquired'
    sys.exit()

# load AsmdMethod_solv_vel_stage.pkl
# ex.: AsmdMethod_vac_02_10.pkl
solvent = my_dir.split('/')[-2].split('.')[1]
vel_dir = my_dir.split('/')[-1]
total_stages = str(int(count)) # '2'
asmd_pkl_name = 'AsmdMethod_%s_%s_%s.pkl' % (solvent,vel_dir,total_stages)
dir_loc_AsmdMethod_pkl = '/'.join(my_dir.split('/')[0:-2])
asmd_pkl = os.path.join(dir_loc_AsmdMethod_pkl,asmd_pkl_name)
sys.path.append(dir_loc_AsmdMethod_pkl)
from asmd.asmdwork import *
c_asmd = pickle.load(open(asmd_pkl,'r'))

print dir(c_asmd)

vel  = c_asmd.v
dist = c_asmd.dist
ts   = c_asmd.ts
path_seg   = c_asmd.path_seg
path_svel  = c_asmd.path_svel
path_vel   = c_asmd.path_vel
path_steps = c_asmd.path_steps
dt         = c_asmd.dt
path_v_aps = c_asmd.pv_aps
domain     = np.cumsum(((path_steps*ts)/1000)*path_v_aps)

print vel
print dist
print ts
print path_seg
print path_svel
print path_vel
print path_steps
print dt
print path_v_aps
print domain

highest_work = []

spos=42.0
kb  =-0.001987
temp=300
beta=1/(kb*temp) # 1/kb*T
quota=10*10
quota = []
spos=42.0
#beta=-0.5961
num =str(len(path_steps)).zfill(2)

dirs = []
for i in range(1,int(num)+1):
    dirs.append(str(i).zfill(2))

#_________________________________________________________________________
def iter_pickle(filename):
    with open(filename) as fp:
        while True:
            try:
                entry = pickle.load(fp)
            except EOFError:
                break
            yield entry

def find_seed(seq,seed):
    for dct in seq:
        if dct['seed']==seed:
            return dct

def residue_index(label):
    return int(re.sub("[^0-9]","",label))
def charac_bond2(trajectory,distance_target):
    acc_count_frames = []
    for frame in trajectory:
        acc_count = 0
        for bond in frame:
            #distance = residue_index(bond[2])-residue_index(bond[3])
            distance = bond[0]-bond[1]
            if distance == distance_target:
                acc_count += 1
        acc_count_frames.append(acc_count)
    return acc_count_frames

def plot_hb(avgB,st,color,lw):
    phase = int(st)-1
    if (st=='01'):
        d = np.linspace(spos,spos+domain[phase],avgB.shape[0])
        ax.plot(d,avgB,color,label="hydrogen bonds",linewidth=lw)
        ax.plot(d[-1],avgB[-1],'.r')
        hbdata=np.transpose([d,avgB])
        np.savetxt('%s_3.dat' % st,hbdata,fmt=['%3.4f','%3.11f'],delimiter=' ')
    elif st !='01':
        d = np.linspace(spos+domain[phase-1],spos+domain[phase],avgB.shape[0])
        ax.plot(d,avgB,color,linewidth=lw)
        ax.plot(d[0],avgB[0],'.r')
        hbdata=np.transpose([d,avgB])
        np.savetxt('%s_3.dat' % st,hbdata,fmt=['%3.4f','%3.11f'],delimiter=' ')


#_________________________________________________________________________

def pack(stage):
    seed_bond={}
    wght_bond={}
    B_list2={}
    W_list2={}
    wrk_pkl={}
    wrk_pkl=pickle.load(open('%s-sfwf.pkl' % stage,'rb'))
    for path in glob(os.path.join(my_dir,'%s/*/*-hb_pr*pr*.pkl.*' % stage)):
        seed = path.split('.')[-1]
        sample_i = pickle.load(open(path,'rb'))
        b_data = np.array(charac_bond2(sample_i, 3))
        seed_bond[seed]=np.array(b_data)
    seeds = wrk_pkl[stage].keys()
    print seeds
    for s in seeds:
        print wrk_pkl[stage][s][1][::,3] # work value for each trajectory/seed
        sample_w = np.exp(wrk_pkl[stage][s][1][::,3]*beta).astype(float)
        sample_b = (seed_bond[s]).astype(float) # bond
        lenf_w = len(sample_w)/100
        lenf_b = len(sample_b)/100
        B_list=[]
        W_list=[]
        H_list=[]
        for b in range(len(sample_b)): # for each frame
            wv = int(((b+1)/lenf_b)*lenf_w)
            #sum_B=(sample_b[b]*np.exp(beta*sample_w[wv]))
            sum_B=sample_b[b]*sample_w[wv] # Bond*exp(-beta*W)
            #sum_W=(np.exp(beta*sample_w[wv]))
            sum_W=sample_w[wv] # exp(-beta*W)
            B_list.append(sum_B) # exponential value for each frame (bond)
            W_list.append(sum_W)  # exponential value for each frame (work)
            H_list.append(sample_b[b])
        #plot_hb_bluedot(np.array(H_list),stage,'b-',0.1)
        B_list2[s]=B_list
        W_list2[s]=W_list
    avg_B=np.sum(B_list2.values(),axis=0)/np.sum(W_list2.values(),axis=0)
    plot_hb(avg_B,stage,'k-',2)

#print b_data
#trj  = str(len(acc))
#trjl = len(acc[0])

#print trj,trj1

dirs = []
for i in range(1,int(num)+1):
    dirs.append(str(i).zfill(2))
[pack(st) for st in sorted(dirs)]

# matplotlib
ax.set_xlabel('end-to-end distance (A)')
ax.set_ylabel('<$N_H(S_P,S_P)$>')
#ax.plot([23,23],[0,10],'r--',lw=1)
#ax.text(12.5,9.5,'(a) Naive ASMD', weight='bold',fontsize=14)
#plt.title('DA | NAMD | ASMD | Charmm22 \n \All Hydrogen Bonds, intra-protein, vac')
#leg=ax.legend()
plt.yticks(np.arange(0,6,2),fontsize=14)
plt.gca().set_ylim(ymin=-0.5)
plt.gca().set_ylim(ymax=4.5)

ax.minorticks_on()
ax.tick_params(axis='both',bottom=True,top=True,left=True,right=True,direction='in',which='both')
#plt.gca().axes.xaxis.set_ticklabels([])
#plt.xticks(np.arange(15,35,5))
axis = plt.gca().yaxis.set_minor_locator(MultipleLocator(1))
axis = plt.gca().xaxis.set_minor_locator(MultipleLocator(1))
plt.subplots_adjust(left=0.11,right=0.99,top=0.97,bottom=0.13)
fig.set_size_inches(7.12,4.4)
plt.legend(['2K6O'])
plt.draw()

texdir = os.path.join(('/'.join(my_dir.split('/')[0:-2])), \
                      'tex_%s/fig_bond' % my_dir.split('/')[-3])
if not os.path.exists(texdir): os.makedirs(texdir)
plotname = 'hbond_pp3_2K6O'
plt.savefig('%s/%s.png' % (texdir,plotname),format='png',dpi=1200)
plt.savefig('%s/%s.eps' % (texdir,plotname),format='eps',dpi=1200)
