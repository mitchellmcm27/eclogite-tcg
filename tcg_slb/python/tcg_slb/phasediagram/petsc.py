import io
import numpy as np
from scipy.integrate import solve_ivp
from contextlib import redirect_stdout, redirect_stderr
import time
import inspect

from ..base import *
from .base import BasePDReactiveODE

from mpi4py import MPI as mpi

PETSC_AVAIL = True
try: 
    from petsc4py import PETSc as petsc
except ImportError:
    PETSC_AVAIL = False

class PETScPDReactiveODE(BasePDReactiveODE):
    class Solution:
        def __init__(self,n):
            self.n = n
            self.i = np.empty(0,dtype=int)
            self.t = np.empty(0,dtype=float)
            self.y = np.empty((n,0),dtype=float)
        def update(self,i,t,y):
            self.i = np.hstack((self.i,i))
            self.t = np.hstack((self.t,t))
            self.y = np.hstack((self.y,y.reshape(self.n,1)))

    comm = None

    ts   = None
    snes = None
    ksp  = None
    pc   = None

    JFMat = None
    JGmat = None
    Fres  = None
    Gres  = None

    print_norms = None

    def __init__(self,rxn):
        super().__init__(rxn)
        
        self.comm = mpi.COMM_WORLD

        self.setup_ts()
        self.setup_mats()
        #self.setup_bounds()

    def setup_ts(self):
        self.ts = petsc.TS().create(comm=self.comm)
        self.ts.setType(self.ts.Type.BDF)
        self.ts.setTime(0.0)
        self.ts.setTimeStep(1.e-6)
        self.ts.setMaxTime(10.0)
        self.ts.setMaxSteps(100000000)
        self.ts.setMaxSNESFailures(-1)
        self.ts.setExactFinalTime(self.ts.ExactFinalTime.MATCHSTEP)
        self.ts.setMonitor(self.TSCustomMonitor)
        
        self.snes = self.ts.getSNES()
        self.snes.setTolerances(max_it=10)
        self.snes.setMonitor(self.SNESCustomMonitor)
        self.snes.setType('vinewtonrsls')
        
        self.ksp = self.snes.getKSP()
        self.ksp.setType(self.ksp.Type.PREONLY)
        
        self.pc = self.ksp.getPC()
        self.pc.setType(self.pc.Type.LU)
        self.pc.setFactorSolverType("mumps")

        self.ts.setFromOptions()
        
    def setup_mats(self):
        self.JFmat = petsc.Mat().create(comm=self.comm)
        sizes = (self.I+self.K, self.I+self.K)
        self.JFmat.setSizes(sizes)
        self.JFmat.setType('aij')
        self.JFmat.setUp()
        
        self.JGmat = petsc.Mat().create(comm=self.comm)
        sizes = (self.I+self.K, self.I+self.K)
        self.JGmat.setSizes(sizes)
        self.JGmat.setType('aij')
        self.JGmat.setUp()
        
        self.Fres = self.JFmat.createVecRight()
        
        self.Gres = self.JGmat.createVecRight()
        
        self.ts.setIFunction(self.FormIFunction, self.Fres)
        self.ts.setIJacobian(self.FormIJacobian, self.JFmat, self.JFmat)
        self.ts.setRHSFunction(self.FormRHSFunction, self.Gres)
        self.ts.setRHSJacobian(self.FormRHSJacobian, self.JGmat, self.JGmat)

    def setup_bounds(self):
        ub = self.Fres.duplicate()
        ub.setArray(np.ones(self.I+self.K))
        lb = self.Fres.duplicate()
        lb.setArray(np.zeros(self.I+self.K))
        self.snes.setVariableBounds(lb, ub)

    def TSCustomMonitor(self,ts,i,t,u):
        dt = ts.getTimeStep()
        print('{} TS dt {} time {}'.format(i,dt,t))
        y = np.array(u[:].tolist())
        self.sol.update(i,t,y)
        
    def SNESCustomMonitor(self,snes,n,norm):
        print('    {} SNES Function norm {}'.format(n,norm))
        
    def solve(self,T,p,mi0,Cik0,end,**kwargs):
        super().set_initial_params(T,p,mi0,Cik0,**kwargs)

        method = kwargs.get('method', 'bdf')
        rtol   = kwargs.get('rtol', 1.e-5)
        atol   = kwargs.get('atol', 1.e-9)
        snes_rtol = kwargs.get('snes_rtol', 1.e-6)
        snes_atol = kwargs.get('snes_atol', 1.e-6)
        snes_maxit = kwargs.get('snes_maxit', 10)
        maxsteps = kwargs.get('maxsteps', 100000000)
        
        self.print_norms = kwargs.get('print_norms', False)

        self.sol = self.Solution(self.I+self.K)

        u0 = np.concatenate((mi0,Cik0))
        x = self.Gres.duplicate()
        x.setArray(u0)
        self.ts.setType(method)
        self.ts.setMaxTime(end)
        self.ts.setTolerances(rtol=rtol,atol=atol)
        self.ts.setMaxSteps(maxsteps)
        self.snes.setTolerances(rtol=snes_rtol, atol=snes_atol, max_it=snes_maxit)
        try:
            so = io.StringIO()
            se = io.StringIO()
            with redirect_stdout(so), redirect_stderr(se):
                tic = time.perf_counter()
                self.ts.solve(x)
                toc = time.perf_counter()
                self.stime = toc-tic
            self.stdout = so.getvalue()
            self.stderr = se.getvalue()
            flag = self.rxn.check_coder_error()
            if flag!=0:
                self.sol = None
                self.excstr = repr(flag)
                self.rxn.reset_coder_error()
        except Exception as e:
            self.sol = None
            self.excstr = repr(e)

    def FormIFunction(self, ts, t, u, udot, f):
        print('  In FormIFunction')
        f.setArray(udot)
        if self.print_norms:
            print('    FormIFunction: 2-norm u = {}'.format(u.norm(norm_type=1)))
            print('    FormIFunction: inf-norm u = {}'.format(u.norm(norm_type=3)))
            print('    FormIFunction: 2-norm u_t = {}'.format(udot.norm(norm_type=1)))
            print('    FormIFunction: inf-norm u_t = {}'.format(udot.norm(norm_type=3)))
            print('    FormIFunction: 2-norm f = {}'.format(f.norm(norm_type=1)))
            print('    FormIFunction: inf-norm f = {}'.format(f.norm(norm_type=3)))
        
    def FormIJacobian(self, ts, t, u, udot, shift, fmat, pfmat):
        print('  In FormIJacobian')
        print("    a (shift) = {}".format(shift))
        x = self.Gres.duplicate()
        x.setArray(shift*np.ones(self.I+self.K))
        fmat.zeroEntries()
        fmat.setDiagonal(x)
        fmat.assemble()
        if self.print_norms:
            print('    FormIJacobian: 2-norm u = {}'.format(u.norm(norm_type=1)))
            print('    FormIJacobian: inf-norm u = {}'.format(u.norm(norm_type=3)))
            print('    FormIJacobian: 2-norm u_t = {}'.format(udot.norm(norm_type=1)))
            print('    FormIJacobian: inf-norm u_t = {}'.format(udot.norm(norm_type=3)))
            print('    FormIJacobian: Frobenius norm A = {}'.format(fmat.norm(norm_type=2)))
            print('    FormIJacobian: inf-norm A = {}'.format(fmat.norm(norm_type=3)))
            print('    FormIJacobian: Frobenius norm B = {}'.format(pfmat.norm(norm_type=2)))
            print('    FormIJacobian: inf-norm B = {}'.format(pfmat.norm(norm_type=3)))
        return True
        
    def FormRHSFunction(self, ts, t, u, g):
        '''
        assemble rhs of the dimensionaless 1-D reactive equations

        Parameters
        ----------

        t: float
            parameterized distance
        u: array
            solution array [ mi, Cik ]

        '''
        print('  In FormRHSFunction')
        du = super().rhs(t, u)
        g.setArray(du)
        if self.print_norms:
            print('    FormRHSFunction: 2-norm u = {}'.format(u.norm(norm_type=1)))
            print('    FormRHSFunction: inf-norm u = {}'.format(u.norm(norm_type=3)))
            print('    FormRHSFunction: 2-norm f = {}'.format(g.norm(norm_type=1)))
            print('    FormRHSFunction: inf-norm f = {}'.format(g.norm(norm_type=3)))


    def FormRHSJacobian(self, ts, t, u, gmat, pgmat):
        '''
        return rhs of the dimensionaless 1-D reactive equations

        Parameters
        ----------

        t: float
            parameterized distance
        u: array
            solution array [ mi, Cik ]

        '''
        print('  In FormRHSJacobian')
        dudu = super().jac(t, u)
        idx = list(range(self.I+self.K))
        gmat.setValues(idx, idx, dudu)
        gmat.assemble()
        if self.print_norms:
            print('    FormRHSJacobian: 2-norm u = {}'.format(u.norm(norm_type=1)))
            print('    FormRHSJacobian: inf-norm u = {}'.format(u.norm(norm_type=3)))
            print('    FormRHSJacobian: Frobenius norm A = {}'.format(gmat.norm(norm_type=2)))
            print('    FormRHSJacobian: inf-norm A = {}'.format(gmat.norm(norm_type=3)))
            print('    FormRHSJacobian: Frobenius norm B = {}'.format(pgmat.norm(norm_type=2)))
            print('    FormRHSJacobian: inf-norm B = {}'.format(pgmat.norm(norm_type=3)))
        return True

if not PETSC_AVAIL:
    class PETScPDReactiveODE(BasePDReactiveODE):
        def __init__(self,rxn):
            raise Exception("PETSc not available")

