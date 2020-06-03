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
parser.add_argument('--nmol', type=int, default=15, help='the `traj` for Local minima')
parser.add_argument('--nstep', type=int, default=30, help='the `traj` for Local minima')
args = parser.parse_args()
######## local minima #########
dr = args.dr
nmol=args.nmol
np.random.seed(2015)

outputdir='h2o_' + str(nmol)
try: 
    os.mkdir(outputdir)
except: 
    print(outputdir, ' FOLDER EXISTS')

for i in range(args.nstep):
    w5=tip4pcluster2(nmol, 5).water()
    ftraj = outputdir + '/lowest_basinm' + str(nmol) +"_"  +str(dr)[2:]  + str(i) +'.traj'

    BH = BasinHoppingm(w5,temperature=300 * kB,
                                     dr=dr,
                                     trajectory=ftraj, logfile= outputdir + '/basin' + str(dr)[2:] + str(i) + '.log')

    BH.run(200)
    Emin, smin = BH.get_minimum()
    print(Emin)