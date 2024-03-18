#!/usr/bin/env python3

import sys, os
import numpy as np
import petsc4py
petsc4py.init(sys.argv)

sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir, 'python'))
from tcg_slb.phasediagram.fenicsx import PDReactiveODE, TSProblem
from tcg_slb.base import *

pv = repr(sys.version_info.major)+'.'+repr(sys.version_info.minor)
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir, 'database', 'install', 'MgFeSiO4_all_slb_rx', 'lib', 'python'+pv, 'site-packages'))

import py_MgFeSiO4_all_slb_rx as tcgdb
rxn = tcgdb.MgFeSiO4_all_slb_rx()
rxn.report()


problem = PDReactiveODE(rxn)
ode = TSProblem(problem)
        
ode.solve(print_norms=False, ts_max_time=1.0)
