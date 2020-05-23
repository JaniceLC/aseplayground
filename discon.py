import os 
import argparse
from ase.io.trajectory import Trajectory


import numpy as np
import pandas as pd
import networkx as nx

from tools.tip4p_cluster import tip4pcluster2
from tools.tip4p_cluster import add_tip4p_const
from tools.tip4p_cluster import prepare_graph


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

# load minima 
traj = args.lm 
empty = os.path.getsize(ftraj) == 0
if not empty:
    traj = io.Trajectory(ftraj, 'r')
    minima = [atoms for atoms in traj]
else:
    print('no minima founded')
    minima = []
nminima = len(minima)
print('# minima: ', nminima)
minima=minima[-args.nmin:]
print('# minima used: ', len(minima))
LM = {}
for i in range(nminima):
    LM[str(i)] = add_tip4p_const(minima[i])

lmdf=pd.DataFrame(columns=['m1', 'E1'])
E = []
for lm in LM.values():
    add_tip4p_const(lm)
    E.append(lm.get_potential_energy())
lmdf['E1']=E
lmdf['m1']=LM.keys()
lmdf['m2']=LM.keys()
lmdf['E2']=E
trandir=args.tr + "/"

# load transition states
transition_states=[name for name in os.listdir(trandir) if os.path.isfile(trandir+name) and os.path.getsize(trandir+name)>0]


g, new= prepare_graph(transition_states, lmdf, trandir)


# plot 
from discon.disconnectivity import DisconnectivityGraph
dg=DisconnectivityGraph(g, csv_map=new, energy_csv=lmdf)
dg.calculate()
import matplotlib.pyplot as plt

dg.plot(title=args.figtitle)
plt.savefig(args.figname+'.png')


