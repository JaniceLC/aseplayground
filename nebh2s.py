import os
import argparse
import math 

from ase import Atoms, io
from ase import units
from math import cos, sin, pi
from ase.calculators.tip4p import TIP4P, epsilon0, sigma0, rOH, angleHOH
from ase.constraints import FixInternals
from ase.neb import NEB
from ase.optimize import MDMin, BFGS,LBFGS,FIRE
from ase.optimize.basin_mul import BasinHoppingm
import pandas as pd
import numpy as np
from tools.h2s_cluster import h2scluster
from tools.h2s_cluster import add_const
from ase.io.trajectory import Trajectory

parser = argparse.ArgumentParser(description='Process input traj files')
parser.add_argument('--lm', type=str, default='h2s_10', help='the `traj` for Local minima')
parser.add_argument('--tr', type=str, default='w10_ts', help='the directory to store ')
parser.add_argument('--nmin', type=int, default=50, help='number of minima used')
args = parser.parse_args()
######## local minima #########
ftraj = args.lm 

from os import listdir
def list_of_files(dir_name, suffix):
    return [f for f in listdir(dir_name) if f.endswith('.' + suffix)]
# load minima 
#ftraj = args.lm 
mydir=args.lm + '/'
minima = list_of_files(args.lm, 'xyz')
nminima = len(minima)
sep = math.floor(nminima/50)
minima=minima[0::sep]
nminima = len(minima)
print('# minima used: ', nminima)

LM = {}
for mini in minima:
    try:
      LM[mini[:-4]]= io.read(mydir+mini)
    except:
      continue
for lm in LM.values():
    add_const(lm)

import itertools
transition=pd.DataFrame(columns=['m1', 'm2', 'Ed', 'images'])
lmlist=LM.keys()
comb = list(itertools.product(lmlist, lmlist))
comb =[list(e) for e in set(frozenset(d) for d in comb)]
dcomb = []
for i in comb:
    if len(i)==2:
        dcomb.append(i)

dcomb=np.array(dcomb)
transition['m1'] = dcomb[:, 0]
#transition['m1']= transition['m1'].apply(lambda x: int(x))
transition['m2'] = dcomb[:, 1]
#transition['m2']= transition['m2'].apply(lambda x: int(x))
transition = transition[transition['m1']!=transition['m2']]
print(transition.head())
npdf=transition.values



########### neb 
import signal
class TimedOutExc(Exception):
    pass

def deadline(timeout, *args):
    def decorate(f):
        def handler(signum, frame):
            raise TimedOutExc()

        def new_f(*args):
            signal.signal(signal.SIGALRM, handler)
            signal.alarm(timeout)
            return f(*args)
            signal.alarm(0)

        new_f.__name__ = f.__name__
        return new_f
    return decorate

@deadline(50)
def optimize(neb="", name=''):
    optimizer = FIRE(neb, trajectory=name+'.traj')
    optimizer.run(fmax=100000)
    optimizer.run(fmax=50000)
    optimizer.run(fmax=10000)
    optimizer.run(fmax=5000)
    optimizer.run(fmax=1000)
    optimizer.run(fmax=500)
    optimizer.run(fmax=100)
    optimizer.run(fmax=50)
    optimizer.run(fmax=10)
    optimizer.run(fmax=0.5)
    optimizer.run(fmax=0.1)
    optimizer.run(fmax=0.09)
    optimizer.run(fmax=0.07)
    optimizer.run(fmax=0.05)

import matplotlib.pyplot as plt
from ase.neb import NEBTools
from ase.io import read
def saddlepoint(min1, min2, name):
    images=[LM[min1]]
    images += [LM[min1].copy() for i in range(5)]
    images += [LM[min2]]
    for image in images:
        image.set_calculator(TIP4P())
        add_const(image)
    neb = NEB(images, climb=True, parallel=True) 
    try:
        neb.interpolate()
    except:
        try:
            print("try idpp interpolate")
            neb.interpolate(method="idpp")
        except:
            print('interploation cannot converge')
            return
    optimize(neb, name)
ts = args.tr
os.mkdir(ts) 
records=0
for line in npdf:
    records+=1
    try:
        line[3]=line[0]+'_'+line[1]
        saddlepoint(line[0], line[1], ts+ "/"+line[3])
    except:
        print(line[0], line[1], "time limit")
        continue
transition=pd.DataFrame(npdf, columns=['m1', 'm2', 'Ed', 'images'])
transition.to_csv(ts+'.csv', index=False)
