
import numpy as np
import pytest
import ase.units as units



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

w5=tip4pcluster2(10, 5).water()

original_positions = 1. * w5.get_positions()
print(original_positions)

hop = MinimaHopping(w5,
                    Ediff0=2.5,
                    #T0=400., 
                   mdmin = 5, 
                    minima_traj = 'w10/lmw10.traj', 
                    logfile='w10/hop.log'
                   )
hop(totalsteps=200)

from ase.optimize.minimahopping import MHPlot

mhplot = MHPlot(logname='w10/hop.log')
mhplot.save_figure('md_w10_200step.png')