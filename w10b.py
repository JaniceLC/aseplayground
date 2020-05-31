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
np.random.seed(2015)
w5=tip4pcluster2(10, 2.5).water()

ftraj = 'lowest_basinm_w10.traj'
nmol=10

for GlobalOptimizer in [BasinHoppingm(w5,
                                     temperature=400 * kB,
                                     dr=0.2,
                                     trajectory=ftraj, logfile='w10basin.log')]:
    if isinstance(GlobalOptimizer, BasinHoppingm):
        GlobalOptimizer.run(300)
        Emin, smin = GlobalOptimizer.get_minimum()
    else:
        GlobalOptimizer(totalsteps=10)
        Emin = w5.get_potential_energy()
        smin =w5
    print("N=", nmol, 'minimal energy found', Emin)
 #         ' global minimum:', E_global[nmol])