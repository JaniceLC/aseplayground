import numpy as np
import pytest
import ase.units as units
from ase import Atoms, io
from ase.constraints import FixInternals
import os
from ase.optimize.minimahopping import MinimaHopping
from ase.optimize.basin_mul import BasinHoppingm


from ase.calculators.lj import LennardJones
#from ase.optimize.basin_mul import BasinHoppingm
from ase.io import read
from ase.units import kB
from ase.calculators.tip4p import TIP4P, epsilon0, sigma0, rOH, angleHOH
from math import cos, sin, pi
from ase.calculators.qmmm import (SimpleQMMM, EIQMMM, LJInteractions,
                                  LJInteractionsGeneral)

from matplotlib import pyplot as plt
from tools.tip4p_cluster import tip4pcluster, tip4pcluster2
from tools.tip4p_cluster import add_tip4p_const
import argparse
parser = argparse.ArgumentParser(description='input optimization para')
parser.add_argument('--dr', type=float, default=0.2, help='the `traj` for Local minima')
args = parser.parse_args()
######## local minima #########
dr = args.dr

np.random.seed(2015)
w5=tip4pcluster2(15, 3.5).water()
nmol=10
try: 
    os.mkdir('h2o15')
except: 
    print('h2o15 FOLDER EXISTS')
ftraj = 'h2o15/lowest_basinm_w15.traj'

BH = BasinHoppingm(w5,temperature=400 * kB,
                                     dr=dr,
                                     trajectory=ftraj, logfile='h2o15/w15basin.log')

BH.run(2000)
Emin, smin = BH.get_minimum()
print(Emin)