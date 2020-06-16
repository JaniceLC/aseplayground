import os 
import argparse
from ase.io.trajectory import Trajectory
import math

import numpy as np
import pandas as pd
import networkx as nx

from tools.h2s_cluster import h2scluster
from tools.h2s_cluster import add_const
from tools.h2s_cluster import prepare_graph
from ase.optimize import BFGS,LBFGS


from ase import Atoms, io
from ase import units
from math import cos, sin, pi
from ase.calculators.tip4p import TIP4P, epsilon0, sigma0, rOH, angleHOH
from ase.constraints import FixInternals
parser = argparse.ArgumentParser(description='Process input traj files')
parser.add_argument('--lm', type=str, default='w10/lmw10.traj', help='the `traj` for Local minima')
parser.add_argument('--tr', type=str, default='w10_ts', help='the directory to store ')
parser.add_argument('--nmin', type=int, default=50, help='number of minima used')
parser.add_argument('--figtitle', type=str, default="TIP4P(H2O)10", help='TITLE OF DISCONNECTIVITY PLOTS e.g. TIP4P(H2O)10 ')
parser.add_argument('--figname', type=str, default="TIP4PW10_DIS", help='figure name')


args = parser.parse_args()
from os import listdir
def list_of_files(dir_name, suffix):
    return [f for f in listdir(dir_name) if f.endswith('.' + suffix)]
# load minima 
#ftraj = args.lm 

mydir=args.lm + '/'
minima = list_of_files(args.lm, 'xyz')
nminima = len(minima)
minima = sorted(minima,key=lambda x: int(os.path.splitext(x)[0]))
#minima.sort(key=lambda f: int(filter(str.isdigit, f)))
sep = math.floor(nminima/args.nmin)
minima=minima[0::sep]
nminima = len(minima)
print('# minima used: ', nminima)

LM = {}

for mini in minima:
    print(mini)
    try:
      LM[mini[:-4]]= io.read(mydir+mini)
    except:
      continue
for lm in LM.values():
    add_const(lm)

for lm in LM.values():
    dyn=LBFGS(lm)
    dyn.run(fmax=0.05)


lmdf=pd.DataFrame(columns=['m1', 'E1'])
E = []
for lm in LM.values():
    add_const(lm)
    E.append(lm.get_potential_energy())
lmdf['E1']=E
lmdf['m1']=LM.keys()
lmdf['m2']=LM.keys()
lmdf['E2']=E
trandir=args.tr + "/"

# load transition states
transition_states=[name for name in os.listdir(trandir) if os.path.isfile(trandir+name) and os.path.getsize(trandir+name)>0]
print(len(transition_states), 'transition paths finded')

g, new= prepare_graph(transition_states, lmdf, trandir)


# plot 
from discon.disconnectivity import DisconnectivityGraph
dg=DisconnectivityGraph(g, csv_map=new, energy_csv=lmdf)
dg.calculate()
import matplotlib.pyplot as plt

dg.plot(title=args.figtitle)
plt.savefig(args.figname+'.png')


