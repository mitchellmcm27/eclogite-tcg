import io
import numpy as np
from scipy.integrate import solve_ivp
from contextlib import redirect_stdout, redirect_stderr
import time
import inspect

from ..base import *
from .petsc import PETScPDReactiveODE

from mpi4py import MPI as mpi

FENICS_AVAIL = True
try: 
    import dolfin as df
except ImportError:
    FENICS_AVAIL = False

class FEniCSPDReactiveODE(PETScPDReactiveODE):

    class TCGExpression(df.UserExpression):
        rxn = None
        shape = None
        func = None
        T_i = None
        p_i = None
        us_i = None
        eps_i = None
        T = None
        p = None
        mi = None
        cik = None
        eps = None
        use_mi = False

        def __init__(self, rxn, name, abbrev, shape, func):

            self.rxn = rxn
            self.I = len(self.rxn.phases())
            self.Kis = [len(self.rxn.phases()[i].endmembers()) for i in range(self.I)]
            self.K = sum(self.Kis)
            self.J = len(self.rxn.nu())
            self.mi = np.zeros(self.I)
            self.cik = self.rxn.zero_C()
            self.svalues = np.zeros(1)
            self.vvalues = np.zeros(self.I+self.K)

            self.shape = shape
            self.func = func
            self.use_mi = len(inspect.getargspec(self.func).args)==5
            super().__init__(degree=0)
            self.rename(name, abbrev)

        def value_shape(self):
            return self.shape

        def eval(self, values, x):
            self.update(x)
            if self.use_mi:
                values[:] = self.func(self.T, self.p, self.cik, self.mi).flatten()
            else:
                values[:] = self.func(self.T, self.p, self.cik).flatten()

        def attach(self, T_i, p_i, us_i, eps_i):
            self.T_i = T_i
            self.p_i = p_i
            self.us_i = us_i
            self.eps_i = eps_i

        def update_T(self, x):
            self.T_i.eval(self.svalues, x)
            self.T = self.svalues[0]
            
        def update_P(self, x):
            self.p_i.eval(self.svalues, x)
            self.p = self.svalues[0]
            
        def update_eps(self, x):
            self.eps_i.eval(self.svalues, x)
            self.eps = self.svalues[0]
            
        def update_mi(self):
            for i in range(self.I):
                self.mi[i] = self.vvalues[i]

        def update_cik(self, x):
            self.update_eps(x)
            sKi = self.I
            for i in range(self.I):
                self.cik[i] = np.clip(self.vvalues[sKi:sKi+self.Kis[i]], self.eps, 1.-self.eps)
                self.cik[i] /= np.sum(self.cik[i])
                sKi += self.Kis[i]

        def update(self, x):
            self.update_T(x)
            self.update_P(x)
            self.us_i.eval(self.vvalues, x)
            self.update_mi()
            self.update_cik(x)

    function_map = {}
    dfunction_map = {}

    def __init__(self,rxn,use_functions=False):
        super(PETScPDReactiveODE,self).__init__(rxn)
        
        self.comm = mpi.COMM_WORLD
        self.use_functions = use_functions

        self.setup_ts()
        self.setup_fe()
        self.setup_mats()
        self.setup_bounds()

    def setup_fe(self):
        self.mesh = df.UnitIntervalMesh(1)
        
        es = [df.FiniteElement("DG", self.mesh.ufl_cell(), 0) for i in range(self.I)]
        es += [df.VectorElement("DG", self.mesh.ufl_cell(), 0, self.Kis[i]) for i in range(self.I)]
        e = df.MixedElement(*es)
        self.V = df.FunctionSpace(self.mesh, e)

        # set these to some dummy values for now, they will be updated later
        self.T_i = df.Constant(666.e6)
        self.T_i.rename("Temperature", "T")
        self.p_i = df.Constant(666.e6)
        self.p_i.rename("Pressure", "p")
        self.Da_i = df.Constant(666.e6)
        self.Da_i.rename("Damkoehler", "Da")
        self.eps_i = df.Constant(666.e6)
        self.eps_i.rename("Epsilon", "eps")
        self.rho0_i = df.Constant(666.e6)
        self.rho0_i.rename("ReferenceDensity", "rho0")
        self.a_i = df.Constant(666.e6)
        self.a_i.rename("Shift", "a")

        self.us_i = df.Function(self.V)
        self.us_i.rename("Solution", "us_i")

        Cik0 = np.ones(self.K)
        for i in range(self.I): Cik0[sum(self.Kis[:i]):sum(self.Kis[:i+1])] /= self.Kis[i]
        self.set_initial_params(1673., 300000.0, np.ones(self.I)/self.I, Cik0)

        us_i_split = df.split(self.us_i)
        self.mi_i = us_i_split[:self.I]
        self.cik_i = us_i_split[self.I:]
        #self.mi_i = [self.us_i.sub(i) for i in range(self.I)]
        #self.cik_i = [self.us_i.sub(i) for i in range(self.I,self.I+self.I)]

        self.us_dot = df.Function(self.V)
        us_dot_split = df.split(self.us_dot)
        self.mi_dot = us_dot_split[:self.I]
        self.cik_dot = us_dot_split[self.I:]
        #self.mi_dot = [self.us_dot.sub(i) for i in range(self.I)]
        #self.cik_dot = [self.us_dot.sub(i) for i in range(self.I,self.I+self.I)]
        
        self.us_t = df.TestFunction(self.V)
        us_t_split = df.split(self.us_t)
        self.mi_t = us_t_split[:self.I]
        self.cik_t = us_t_split[self.I:]

        self.us_a = df.TrialFunction(self.V)
        us_a_split = df.split(self.us_a)
        self.mi_a = us_a_split[:self.I]
        self.cik_a = us_a_split[self.I:]

        rhoi_i = self.TCGExpression(self.rxn, "PhaseDensities", "rhoi", (self.I,), self.rhoi)
        rhoi_i.attach(self.T_i, self.p_i, self.us_i, self.eps_i)
        drhoidus_i = self.TCGExpression(self.rxn, "PhaseDensityDerivatives", "drhoidus", (self.I,self.I+self.K), self.drhoidus)
        drhoidus_i.attach(self.T_i, self.p_i, self.us_i, self.eps_i)
        Gammai_i = self.TCGExpression(self.rxn, "PhaseSources", "Gammai", (self.I,), self.Gammai)
        Gammai_i.attach(self.T_i, self.p_i, self.us_i, self.eps_i)
        dGammaidus_i = self.TCGExpression(self.rxn, "PhaseSourceDerivatives", "dGammaidus", (self.I,self.I+self.K), self.dGammaidus)
        dGammaidus_i.attach(self.T_i, self.p_i, self.us_i, self.eps_i)
        Gammaik_i = self.TCGExpression(self.rxn, "ComponentSources", "Gammaik", (self.K,), self.Gammaik)
        Gammaik_i.attach(self.T_i, self.p_i, self.us_i, self.eps_i)
        dGammaikdus_i = self.TCGExpression(self.rxn, "ComponentSourceDerivatives", "dGammaikdus", (self.K,self.I+self.K), self.dGammaikdus)
        dGammaikdus_i.attach(self.T_i, self.p_i, self.us_i, self.eps_i)

        if self.use_functions:
          e_rhoi_i = df.VectorElement("DG", self.mesh.ufl_cell(), 0, self.I)
          V_rhoi_i = df.FunctionSpace(self.mesh, e_rhoi_i)
          self.rhoi_i = df.Function(V_rhoi_i)
          self.rhoi_i.rename("PhaseDensitiesFunction", "rhoi_i")
          self.function_map[self.rhoi_i] = rhoi_i

          e_drhoidus_i = df.TensorElement("DG", self.mesh.ufl_cell(), 0, (self.I,self.I+self.K))
          V_drhoidus_i = df.FunctionSpace(self.mesh, e_drhoidus_i)
          self.drhoidus_i = df.Function(V_drhoidus_i)
          self.drhoidus_i.rename("PhaseDensityDerivativesFunction", "drhoidus_i")
          self.dfunction_map[self.drhoidus_i] = drhoidus_i

          e_Gammai_i = df.VectorElement("DG", self.mesh.ufl_cell(), 0, self.I)
          V_Gammai_i = df.FunctionSpace(self.mesh, e_Gammai_i)
          self.Gammai_i = df.Function(V_Gammai_i)
          self.Gammai_i.rename("PhaseSourcesFunction", "Gammai_i")
          self.function_map[self.Gammai_i] = Gammai_i

          e_dGammaidus_i = df.TensorElement("DG", self.mesh.ufl_cell(), 0, (self.I,self.I+self.K))
          V_dGammaidus_i = df.FunctionSpace(self.mesh, e_dGammaidus_i)
          self.dGammaidus_i = df.Function(V_dGammaidus_i)
          self.dGammaidus_i.rename("PhaseSourceDerivativesFunction", "dGammaidus_i")
          self.dfunction_map[self.dGammaidus_i] = dGammaidus_i

          e_Gammaik_i = df.VectorElement("DG", self.mesh.ufl_cell(), 0, self.K)
          V_Gammaik_i = df.FunctionSpace(self.mesh, e_Gammaik_i)
          self.Gammaik_i = df.Function(V_Gammaik_i)
          self.Gammaik_i.rename("ComponentSourcesFunction", "Gammaik_i")
          self.function_map[self.Gammaik_i] = Gammaik_i

          e_dGammaikdus_i = df.TensorElement("DG", self.mesh.ufl_cell(), 0, (self.K,self.I+self.K))
          V_dGammaikdus_i = df.FunctionSpace(self.mesh, e_dGammaikdus_i)
          self.dGammaikdus_i = df.Function(V_dGammaikdus_i)
          self.dGammaikdus_i.rename("ComponentSourceDerivativesFunction", "dGammaikdus_i")
          self.dfunction_map[self.dGammaikdus_i] = dGammaikdus_i
        else:
          self.rhoi_i        = rhoi_i
          self.drhoidus_i    = drhoidus_i
          self.Gammai_i      = Gammai_i
          self.dGammaidus_i  = dGammaidus_i
          self.Gammaik_i      = Gammaik_i
          self.dGammaikdus_i = dGammaikdus_i   
     
        v_i = sum([self.mi_i[i]/self.rhoi_i[i] for i in range(self.I)])
        self.F = 0
        self.G = 0
        sKi = 0
        for i in range(self.I):
            self.F += self.mi_t[i]*self.mi_dot[i]*df.dx
            self.G += self.mi_t[i]*self.Da_i*self.rho0_i*self.Gammai_i[i]*v_i*df.dx
            Ki = self.cik_i[i].ufl_shape[0]
            for k in range(Ki):
                GikcGi_i = self.Gammaik_i[sKi+k] - self.cik_i[i][k]*self.Gammai_i[i]
                self.F += self.cik_t[i][k]*self.cik_dot[i][k]*df.dx
                self.G += self.cik_t[i][k]*self.Da_i*self.rho0_i*GikcGi_i*v_i/(self.mi_i[i] + self.eps_i)*df.dx
            sKi += Ki

        cd = {
              self.Gammai_i : self.dGammaidus_i, \
              self.Gammaik_i : self.dGammaikdus_i, \
              self.rhoi_i : self.drhoidus_i, \
             }
        self.JF = self.a_i*df.derivative(self.F, self.us_dot, self.us_a, cd) + df.derivative(self.F, self.us_i, self.us_a, cd)
        self.JG = df.derivative(self.G, self.us_i, self.us_a, cd)

    def setup_mats(self):
        self.JFmat = df.PETScMatrix()
        self.Fres = df.PETScVector()
        Fsysassembler = df.SystemAssembler(self.JF, self.F)
        Fsysassembler.keep_diagonal = True
        Fsysassembler.assemble(self.JFmat, self.Fres)

        self.JGmat = df.PETScMatrix()
        self.Gres = df.PETScVector()
        Gsysassembler = df.SystemAssembler(self.JG, self.G)
        Gsysassembler.keep_diagonal = True
        Gsysassembler.assemble(self.JGmat, self.Gres)

        #self.rxn.reset_coder_error()

        self.ts.setIFunction(self.FormIFunction, self.Fres.vec())
        self.ts.setIJacobian(self.FormIJacobian, self.JFmat.mat(), self.JFmat.mat())
        self.ts.setRHSFunction(self.FormRHSFunction, self.Gres.vec())
        self.ts.setRHSJacobian(self.FormRHSJacobian, self.JGmat.mat(), self.JGmat.mat())

    def setup_bounds(self):
        ub = df.PETScVector(self.comm, self.us_i.vector().size())
        ub.set_local(np.ones(ub.local_size()))
        ub.apply("insert")
        lb = df.PETScVector(self.comm, self.us_i.vector().size())
        lb.zero()
        self.snes.setVariableBounds(lb.vec(), ub.vec())

    def set_initial_params(self,T,p,mi0,Cik0,**kwargs):
        super().set_initial_params(T,p,mi0,Cik0,**kwargs)
        self.T_i.assign(self.T)
        self.p_i.assign(self.p)
        self.Da_i.assign(self.Da)
        self.eps_i.assign(self.eps)
        self.rho0_i.assign(self.rho0)
        u0 = df.Constant(np.concatenate((self.mi0,self.Cik0)))
        self.us_i.interpolate(u0)

    def update_functions(self):
        for f,e in self.function_map.items():
            f.interpolate(e)

    def update_dfunctions(self):
        for f,e in self.dfunction_map.items():
            f.interpolate(e)

    def solve(self,T,p,mi0,Cik0,end,**kwargs):
        self.set_initial_params(T,p,mi0,Cik0,**kwargs)

        method = kwargs.get('method', 'bdf')
        rtol   = kwargs.get('rtol', 1.e-5)
        atol   = kwargs.get('atol', 1.e-9)
        snes_rtol = kwargs.get('snes_rtol', 1.e-6)
        snes_atol = kwargs.get('snes_atol', 1.e-6)
        snes_maxit = kwargs.get('snes_maxit', 10)
        maxsteps = kwargs.get('maxsteps', 100000000)
        
        self.print_norms = kwargs.get('print_norms', False)

        self.sol = self.Solution(self.I+self.K)

        u0 = df.Constant(np.concatenate((mi0,Cik0)))
        self.us_i.interpolate(u0)
        x = df.PETScVector(self.us_i.vector().vec())
        self.ts.setType(method)
        self.ts.setMaxTime(end)
        self.ts.setTolerances(rtol=rtol,atol=atol)
        self.ts.setMaxSteps(maxsteps)
        self.snes.setTolerances(rtol=snes_rtol, atol=snes_atol, max_it=snes_maxit)
        #try:
        so = io.StringIO()
        se = io.StringIO()
        #with redirect_stdout(so), redirect_stderr(se):
        tic = time.perf_counter()
        self.ts.solve(x.vec())
        reason = self.ts.getConvergedReason()
        print("  TSConvergedReason: ", reason)
        toc = time.perf_counter()
        self.stime = toc-tic
        #
        self.stdout = so.getvalue()
        self.stderr = se.getvalue()
        flag = self.rxn.check_coder_error()
        #if flag!=0:
        #    self.sol = None
        #    self.excstr = repr(flag)
        #    self.rxn.reset_coder_error()
        #except Exception as e:
        #    self.sol = None
        #    self.excstr = repr(e)

    def TSCustomMonitor(self,ts,i,t,u):
        dt = ts.getTimeStep()
        print('{} TS dt {} time {}'.format(i,dt,t))
        uv = df.PETScVector(u)
        self.us_i.vector()[:] = uv
        y = np.zeros(self.I+self.K)
        self.us_i.eval(y, 0.5*np.ones(1))
        self.sol.update(i,t,y)
        
    def FormIFunction(self, ts, t, u, udot, f):
        print('  In FormIFunction')
        udotv = df.PETScVector(udot)
        uv = df.PETScVector(u)
        
        self.us_dot.vector()[:] = udotv
        self.us_i.vector()[:] = uv
        
        rhs = df.PETScVector(f)
        sysassembler = df.SystemAssembler(self.JF, self.F)
        sysassembler.assemble(rhs)

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

        self.a_i.assign(shift)

        udotv = df.PETScVector(udot)
        uv = df.PETScVector(u)
        
        self.us_dot.vector()[:] = udotv
        self.us_i.vector()[:] = uv

        matrix = df.PETScMatrix(fmat)
        sysassembler = df.SystemAssembler(self.JF, self.F)
        sysassembler.keep_diagonal = True
        sysassembler.assemble(matrix)
        
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
        uv = df.PETScVector(u)
        self.us_i.vector()[:] = uv
        
        self.update_functions()

        rhs = df.PETScVector(g)
        sysassembler = df.SystemAssembler(self.JF, self.G)
        sysassembler.assemble(rhs)

        if self.print_norms:
            print('    FormRHSFunction: 2-norm u = {}'.format(u.norm(norm_type=1)))
            print('    FormRHSFunction: inf-norm u = {}'.format(u.norm(norm_type=3)))
            print('    FormRHSFunction: 2-norm g = {}'.format(g.norm(norm_type=1)))
            print('    FormRHSFunction: inf-norm g = {}'.format(g.norm(norm_type=3)))


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
        uv = df.PETScVector(u)
        self.us_i.vector()[:] = uv

        self.update_functions()
        self.update_dfunctions()

        matrix = df.PETScMatrix(gmat)
        sysassembler = df.SystemAssembler(self.JG, self.G)
        sysassembler.keep_diagonal = True
        sysassembler.assemble(matrix)
        
        if self.print_norms:
            print('    FormRHSJacobian: 2-norm u = {}'.format(u.norm(norm_type=1)))
            print('    FormRHSJacobian: inf-norm u = {}'.format(u.norm(norm_type=3)))
            print('    FormRHSJacobian: Frobenius norm A = {}'.format(gmat.norm(norm_type=2)))
            print('    FormRHSJacobian: inf-norm A = {}'.format(gmat.norm(norm_type=3)))
            print('    FormRHSJacobian: Frobenius norm B = {}'.format(pgmat.norm(norm_type=2)))
            print('    FormRHSJacobian: inf-norm B = {}'.format(pgmat.norm(norm_type=3)))
        return True

if not FENICS_AVAIL:
    class FEniCSPDReactiveODE(PETScPDReactiveODE):
        def __init__(self,rxn):
            raise Exception("FEniCS not available")


