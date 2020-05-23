
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
from ase.calculators.tip4p import TIP4P, epsilon0, sigma0, rOH, angleHOH
from math import cos, sin, pi
from ase.calculators.qmmm import (SimpleQMMM, EIQMMM, LJInteractions,
                                  LJInteractionsGeneral)

from matplotlib import pyplot as plt
from tools.tip4p_cluster import tip4pcluster, tip4pcluster2
from tools.tip4p_cluster import add_tip4p_const
np.random.seed(86)
w5=tip4pcluster2(10, 3).water()

original_positions = 1. * w5.get_positions()
print(original_positions)

hop = MinimaHopping(w5,
                    Ediff0=0.5,
                    #T0=400., 
                   mdmin = 5, 
                   optimizer=GPMin,
                    minima_traj = 'w10/lmw10_compact.traj', 
                    logfile='w10/hop_compact.log'
                   )
hop(totalsteps=500)

from ase.optimize.minimahopping import MHPlot

mhplot = MHPlot(logname='w10/hop_compact.log')
mhplot.save_figure('md_w10_500step_compact.png')