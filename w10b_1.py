import numpy as np
import pytest
import ase.units as units
from ase import Atoms, io
from ase.constraints import FixInternals

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
import os
np.random.seed(2015)
w5=tip4pcluster2(10, 2.5).water()
try: 
    os.mkdir('h2o10')
except: 
    print('h2o10 FOLDER EXISTS')

ftraj = 'h2o10/lowest_basinm_w10_small.traj'
nmol=10

BH = BasinHoppingm(w5,temperature=400 * kB,
                                     dr=0.1,
                                     trajectory=ftraj, logfile='h2o10/w10basin_small.log')

BH.run(800)
Emin, smin = BH.get_minimum()
print(Emin)