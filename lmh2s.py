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
np.random.seed(2029)
w5=h2scluster(5, 3).water()

original_positions = 1. * w5.get_positions()
print(original_positions)
try: 
    os.mkdir('h2s5')
except: 
    os.rmdir('h2s5')
    os.mkdir('h2s5')
hop = MinimaHopping(w5,
                    Ediff0=0.5,
                    #T0=400., 
                   mdmin = 5, 
                   minima_traj = 'h2s5/h2s5.traj', 
                   logfile='h2s5/h2s5.log'
                   )
hop(totalsteps=20)

from ase.optimize.minimahopping import MHPlot

mhplot = MHPlot(logname='h2s5/h2s5.log')
mhplot.save_figure('h2s5/h2s5.png')