import os
import numpy as np
import ase.units as units
from ase.optimize import BFGS, GPMin
from ase import Atoms, io
from ase.constraints import FixInternals

from ase.optimize.minimahopping import MinimaHopping

from ase.calculators.lj import LennardJones
#from ase.optimize.basin_mul import BasinHoppingm
from ase.io import read
from ase.units import kB
from ase.calculators.h2s import H2S, epsilon0, sigma0, rSH, angleHSH
from math import cos, sin, pi
from ase.calculators.qmmm import (SimpleQMMM, EIQMMM, LJInteractions,
                                  LJInteractionsGeneral)

from matplotlib import pyplot as plt
from tools.h2s_cluster import h2scluster
from tools.h2s_cluster import add_const
np.random.seed(2030)
w5=h2scluster(5, 3).water()

original_positions = 1. * w5.get_positions()
print(original_positions)
print('try local minimization')
calc=H2S()
from ase.optimize import BFGS, LBFGS
from ase.optimize.basin_mul import BasinHoppingm
opt = BFGS(w5, maxstep=0.02,
               trajectory='w5' + calc.name + '.traj', 
           logfile='w5' + calc.name + 'd.log')

opt.run(50000)
opt.run(5000)
opt.run(500)
opt.run(100)
opt.run(10)
opt.run(1)

print('local min finished')

try: 
    os.mkdir('h2s5')
except: 
    print('folder h2s5 exists')

BH = BasinHoppingm(w5,
                                     temperature=400 * kB,
                                     dr=0.2,
                                     trajectory='h2s5/h2s5.traj', 
                                     logfile='h2s5/h2s5.log')

BH.run(800)
Emin, smin = BH.get_minimum()
print(Emin)
