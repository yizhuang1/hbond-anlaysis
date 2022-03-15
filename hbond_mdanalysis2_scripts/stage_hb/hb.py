#!/usr/bin/env python
import MDAnalysis
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
from sys import argv
import numpy as np
import os,sys,pickle
import pandas as pd
import matplotlib.pyplot as plt

#______________universe________________________________________________________
u = MDAnalysis.Universe('../../../../00.struc/03.exp/00.psf','daOut.dcd',\
                        permissive=True)

def analyze_bond(univ,seg1,seg2):
#    try:
        name1=seg1.replace(' ','')
        name2=seg2.replace(' ','')
        hbonds = HBA(universe=u,between=[seg1, seg2], d_a_cutoff=4.0, d_h_a_angle_cutoff=140.0, update_selections=True)
        protein_hydrogens_sel = hbonds.guess_hydrogens("protein")
        protein_acceptors_sel = hbonds.guess_acceptors("protein")

        water_hydrogens_sel = "resname TIP3 and name H1 H2"
        water_acceptors_sel = "resname TIP3 and name OH2"

        hbonds.hydrogens_sel = f"({protein_hydrogens_sel}) or ({water_hydrogens_sel} and around 10 not resname TIP3)"
        hbonds.acceptors_sel = f"({protein_acceptors_sel}) or ({water_acceptors_sel} and around 10 not resname TIP3)"
        hbonds.run()
        FRAME = 0
        DONOR = 1
        HYDROGEN = 2
        ACCEPTOR = 3
        DISTANCE = 4
        ANGLE = 5
        df = pd.DataFrame(hbonds.results.hbonds[:, :DISTANCE].astype(int),\
                  columns=["Frame",\
                           "Donor_ix",\
                           "Hydrogen_ix",\
                           "Acceptor_ix",])

        df["Distances"] = hbonds.results.hbonds[:, DISTANCE]
        df["Angles"] = hbonds.results.hbonds[:, ANGLE]
        df["Donor resname"] = u.atoms[df.Donor_ix].resnames
        df["Acceptor resname"] = u.atoms[df.Acceptor_ix].resnames
        df["Donor resid"] = u.atoms[df.Donor_ix].resids
        df["Acceptor resid"] = u.atoms[df.Acceptor_ix].resids
        df["Donor name"] = u.atoms[df.Donor_ix].names
        df["Acceptor name"] = u.atoms[df.Acceptor_ix].names
        df.to_csv('%s-hb_%s_%s.csv.%s' %(sys.argv[1],name1,name2,sys.argv[2]),index=False)
        dump_data=[]
        for n in range(len(u.trajectory)):
            df2=df.loc[df['Frame']==n]
            sublist=[]
            for index, row in df2.iterrows():
                sublist.append([row['Donor resid'],row['Acceptor resid'],str(row['Donor resname']+':'+row['Donor name']).replace("u'", ""),str(row['Acceptor resname']+":"+row['Acceptor name']).replace("u'", ""),row['Distances'],row['Angles']])
            dump_data.append(sublist)
        print(len(dump_data))
        pickle.dump(dump_data,open('%s-hb_%s_%s.pkl.%s' %(sys.argv[1],name1,name2,sys.argv[2]),'wb'),protocol=2)
        np.savetxt('%s-hb_%s_%s.dat.%s' %(sys.argv[1],name1,name2,sys.argv[2]),np.transpose([hbonds.times,hbonds.count_by_time()]),fmt=['%.2f', '%.2f'],delimiter=' ')
#    except:
#        pass

#__analyze__bonds______________________________________________________________
analyze_bond(u,'protein','protein')
analyze_bond(u,'protein','resname TIP3')
#analyze_bond(u,'protein','segid WT1')
