#!/usr/bin/env python3

import sys, os
import numpy as np
import petsc4py
petsc4py.init(sys.argv)

sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir, 'python'))
from tcg_slb.phasediagram.fenics import FEniCSPDReactiveODE
from tcg_slb.base import *

pv = repr(sys.version_info.major)+'.'+repr(sys.version_info.minor)
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir, 'database', 'install', 'MgFeSiO4_all_slb_rx', 'lib', 'python'+pv, 'site-packages'))

import py_MgFeSiO4_all_slb_rx as tcgdb
rxn = tcgdb.MgFeSiO4_all_slb_rx()
rxn.report()

odef = FEniCSPDReactiveODE(rxn)

# initial temperature, pressure and phase volume fraction
Ti = 1673.                # Kelvin
pi = 150000.         # bars
Ci0 = [0.9, 0.1]
i0 =  0                   # initial phase index
end = 1

Fi = np.zeros(odef.I)
Fi[i0] = 1.
Cik = np.zeros(odef.K)
for i in range(odef.I):
    if odef.Kis[i] == 1:
        Cik[sum(odef.Kis[:i]):sum(odef.Kis[:i+1])] = 1.
    else:
        Cik[sum(odef.Kis[:i]):sum(odef.Kis[:i+1])][:2] = np.asarray(Ci0)
        
odef.solve(Ti,pi,Fi,Cik,end,print_norms=True)
