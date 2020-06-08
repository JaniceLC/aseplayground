
from ase import io
from math import cos, sin, pi
from ase import Atoms, io
import numpy as np
from ase.constraints import FixInternals
# from ase.calculators.tip4p import TIP4P, epsilon0, sigma0, rOH, angleHSH
from ase.calculators.h2s import H2S, epsilon0, sigma0, rSH, angleHSH

import ase.units as units
import numpy as np
import matplotlib.pyplot as plt
from ase.neb import NEBTools
from ase.io import read
import networkx as nx
import pandas as pd

class h2scluster:
    def __init__(self, nmol=5, size=3):
        r = rSH
        a = angleHSH * pi / 180
        
        # From https://doi.org/10.1063/1.445869
        eexp = 6.24 * units.kcal / units.mol
        dexp = 2.75
        aexp = 46
        D = np.linspace(2.5, 3.5, 30)
        x = angleHSH * np.pi / 180 / 2
  
        pos=  [[0.00000000  ,  0.00000000,    0.10357400],
                [0.00000000  ,  0.97644400 ,  -0.82858900],
                [0.00000000  , -0.97644400 ,  -0.82858900]]
        pos = pos*nmol
        calc = H2S()
        self.h2o = Atoms('SH2'*nmol, pos)
        disp = np.random.uniform(-size, size, (nmol, 3))
        final_dis=[]
        for i in range(nmol):
            d=disp[i].tolist()
            for i in range(3):
                final_dis.append(d)

        self.h2o.set_positions(pos +np.array(final_dis ))
        #vol = ((18.01528 / 6.022140857e23) / (0.9982 / 1e24))**(1 / 3.)
        #self.h2o.set_cell((1.1*r , 1.1*r , 1.1*r ))
        #nmol = int(self.h2o.positions.shape[0] / 3)
        bonds = [(3 * i, 3 * i + 1) for i in range(nmol)]
        for i in range(nmol):
            bonds.append((3 * i, 3 * i + 2))
        bonds_rSH = []
        for bond in bonds:
            bonds_rSH.append((rSH, (bond)))
        angles = [(a, (3 * i + 2, 3 * i, 3 * i + 1)) for i in range(nmol)]
        self.h2o.constraints = FixInternals(
            bonds=bonds_rSH,
            angles=angles, epsilon=epsilon0)
        calc = H2S()
        self.h2o.calc = calc
    def water(self):
        return self.h2o



def add_const(water):
    nmol = int(water.positions.shape[0] / 3)
    a = angleHSH * pi / 180
    bonds = [(3 * i, 3 * i + 1) for i in range(nmol)]
    for i in range(nmol):
        bonds.append((3 * i, 3 * i + 2))
    bonds_rSH = []
    for bond in bonds:
        bonds_rSH.append((rSH, (bond)))
    angles = [(a, (3 * i + 2, 3 * i, 3 * i + 1)) for i in range(nmol)]
    water.constraints = FixInternals(
        bonds=bonds_rSH,
        angles=angles, epsilon=epsilon0)
    calc = H2S()
    water.calc = calc
    return water 


def prepare_graph(transition_states, lmdf, trandir, log='@-7:'):
    Ed=[]
    for ts in transition_states:
        outputs = read(trandir + ts + log)
        for output in outputs:
            add_const(output)
        nebtools = NEBTools(outputs)
        Ef, dE = nebtools.get_barrier()
        Ed.append([ts[:-5], Ef, Ed])
    Ed=pd.DataFrame(Ed, columns=['images', 'Ed', 'deltaE'])
    mins = Ed['images'].apply(lambda x: x.split("_"))
    min1=[]
    min2=[]
    revers=[]
    for line in mins:
        min1.append(line[0])
        min2.append(line[1])
        revers.append(line[1]+"_"+line[0])
    db = pd.DataFrame(columns=["m1", "m2", "Ed", "images"])


    db['m1'] = min1
    db['m2'] = min2
    db['Ed'] = Ed['Ed']
    db['images'] = Ed['images']
    db = db.merge(lmdf[['m2', 'E2']], how='left', on='m2')


    db = db.merge(lmdf[['m1', 'E1']], how='left', on='m1')
    db['ET'] = db['Ed'] + db['E1']
    db2 = pd.DataFrame(
        columns=["m1", "m2", "Ed", "images", "E2", "E1", "ET"])


    db2['m1'] = db["m2"]
    db2['m2'] = db["m1"]
    db2['ET'] = db['ET']
    db2['Ed'] = db['Ed']
    db2['E1'] = db["E2"]
    db2['E2'] = db["E1"]
    db2['images'] = revers
    new = pd.concat([db, db2])


    new = new[new["ET"] <= 0]
    new = new[new["E1"] <= new["ET"]]
    new = new[new["E2"] <= new["ET"]]
    dbnp = new.values


    g = nx.Graph()
    g.add_nodes_from(new['m1'])
    g.add_nodes_from(new['m2'])
    for line in dbnp:
        line[3] = line[0][4:] + '_' + line[1][4:]
        if line[2] != np.nan:
            g.add_edge(line[0], line[1], ts=[line[0], line[1], line[6]])
    return g, new
            

