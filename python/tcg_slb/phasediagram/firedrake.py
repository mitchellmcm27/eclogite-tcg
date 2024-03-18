import io
import numpy as np
from contextlib import redirect_stdout, redirect_stderr
import time
import inspect

from ..base import *
from .petsc import PETScPDReactiveODE

from mpi4py import MPI as mpi

FIREDRAKE_AVAIL = True
try: 
    import firedrake as fd
except ImportError:
    FIREDRAKE_AVAIL = False

class FiredrakePDReactiveODE(PETScPDReactiveODE):

    class TCGExpression:
        rxn = None
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

        def __init__(self, rxn, func):

            self.rxn = rxn
            self.I = len(self.rxn.phases())
            self.Kis = [len(self.rxn.phases()[i].endmembers()) for i in range(self.I)]
            self.K = sum(self.Kis)
            self.J = len(self.rxn.nu())
            self.mi = np.zeros(self.I)
            self.cik = self.rxn.zero_C()
            self.usv = np.zeros(self.I+self.K)

            self.func = func
            self.use_mi = len(inspect.getargspec(self.func).args)==5

        def evaluate(self, x):
            self.update(x)
            if self.use_mi:
                return self.func(self.T, self.p, self.cik, self.mi).flatten()
            else:
                return self.func(self.T, self.p, self.cik).flatten()

        def attach(self, T_i, p_i, us_i, eps_i):
            self.T_i = T_i
            self.p_i = p_i
            self.us_i = us_i
            self.eps_i = eps_i

        def update_T(self, x):
            self.T = float(self.T_i)
            
        def update_P(self, x):
            self.p = float(self.p_i)
            
        def update_eps(self, x):
            self.eps = float(self.eps_i)
            
        def update_mi(self):
            for i in range(self.I):
                self.mi[i] = self.usv[i]

        def update_cik(self, x):
            self.update_eps(x)
            sKi = self.I
            for i in range(self.I):
                self.cik[i] = np.clip(self.usv[sKi:sKi+self.Kis[i]], self.eps, 1.-self.eps)
                self.cik[i] /= np.sum(self.cik[i])
                sKi += self.Kis[i]

        def update(self, x):
            self.update_T(x)
            self.update_P(x)
            self.usv = self.us_i.vector().array()
            self.update_mi()
            self.update_cik(x)

    function_map = {}
    dfunction_map = {}

    def __init__(self,rxn):
        super(PETScPDReactiveODE,self).__init__(rxn)
        
        self.comm = mpi.COMM_WORLD
        self.use_functions = True # has to be true for firedrake I think

        self.setup_ts()
        self.setup_fe()
        self.setup_mats()
        #self.setup_bounds()

    def setup_ts(self):
        super().setup_ts()

    def setup_fe(self):
        self.mesh = fd.UnitIntervalMesh(1)
        
        Vs = [fd.FunctionSpace(self.mesh, "DG", 0) for i in range(self.I)]
        Vs += [fd.VectorFunctionSpace(self.mesh, "DG", 0, self.Kis[i]) for i in range(self.I)]
        self.V = fd.MixedFunctionSpace(Vs)

        # set these to some dummy values for now, they will be updated later
        self.T_i = fd.Constant(666.e6)
        self.p_i = fd.Constant(666.e6)
        self.Da_i = fd.Constant(666.e6)
        self.eps_i = fd.Constant(666.e6)
        self.rho0_i = fd.Constant(666.e6)
        self.a_i = fd.Constant(666.e6)

        self.us_i = fd.Function(self.V)
        self.us_i.rename("Solution", "us_i")

        Cik0 = np.ones(self.K)
        for i in range(self.I): Cik0[sum(self.Kis[:i]):sum(self.Kis[:i+1])] /= self.Kis[i]
        self.set_initial_params(1673., 300000.0, np.ones(self.I)/self.I, Cik0)

        us_i_split = fd.split(self.us_i)
        self.mi_i = us_i_split[:self.I]
        self.cik_i = us_i_split[self.I:]
        #self.mi_i = [self.us_i.sub(i) for i in range(self.I)]
        #self.cik_i = [self.us_i.sub(i) for i in range(self.I,self.I+self.I)]

        self.us_dot = fd.Function(self.V)
        us_dot_split = fd.split(self.us_dot)
        self.mi_dot = us_dot_split[:self.I]
        self.cik_dot = us_dot_split[self.I:]
        #self.mi_dot = [self.us_dot.sub(i) for i in range(self.I)]
        #self.cik_dot = [self.us_dot.sub(i) for i in range(self.I,self.I+self.I)]
        
        self.us_t = fd.TestFunction(self.V)
        us_t_split = fd.split(self.us_t)
        self.mi_t = us_t_split[:self.I]
        self.cik_t = us_t_split[self.I:]

        self.us_a = fd.TrialFunction(self.V)
        us_a_split = fd.split(self.us_a)
        self.mi_a = us_a_split[:self.I]
        self.cik_a = us_a_split[self.I:]

        self.e_rhoi_i = self.TCGExpression(self.rxn, self.rhoi)
        self.e_rhoi_i.attach(self.T_i, self.p_i, self.us_i, self.eps_i)
        self.e_drhoidus_i = self.TCGExpression(self.rxn, self.drhoidus)
        self.e_drhoidus_i.attach(self.T_i, self.p_i, self.us_i, self.eps_i)
        self.e_Gammai_i = self.TCGExpression(self.rxn, self.Gammai)
        self.e_Gammai_i.attach(self.T_i, self.p_i, self.us_i, self.eps_i)
        self.e_dGammaidus_i = self.TCGExpression(self.rxn, self.dGammaidus)
        self.e_dGammaidus_i.attach(self.T_i, self.p_i, self.us_i, self.eps_i)
        self.e_Gammaik_i = self.TCGExpression(self.rxn, self.Gammaik)
        self.e_Gammaik_i.attach(self.T_i, self.p_i, self.us_i, self.eps_i)
        self.e_dGammaikdus_i = self.TCGExpression(self.rxn, self.dGammaikdus)
        self.e_dGammaikdus_i.attach(self.T_i, self.p_i, self.us_i, self.eps_i)

        if self.use_functions:
          V_rhoi_i = fd.VectorFunctionSpace(self.mesh, "DG", 0, self.I)
          self.rhoi_i = fd.Function(V_rhoi_i)
          self.rhoi_i.rename("PhaseDensitiesFunction", "rhoi_i")
          self.function_map[self.rhoi_i] = self.e_rhoi_i

          V_drhoidus_i = fd.TensorFunctionSpace(self.mesh, "DG", 0, (self.I,self.I+self.K))
          self.drhoidus_i = fd.Function(V_drhoidus_i)
          self.drhoidus_i.rename("PhaseDensityDerivativesFunction", "drhoidus_i")
          self.dfunction_map[self.drhoidus_i] = self.e_drhoidus_i

          V_Gammai_i = fd.VectorFunctionSpace(self.mesh, "DG", 0, self.I)
          self.Gammai_i = fd.Function(V_Gammai_i)
          self.Gammai_i.rename("PhaseSourcesFunction", "Gammai_i")
          self.function_map[self.Gammai_i] = self.e_Gammai_i

          V_dGammaidus_i = fd.TensorFunctionSpace(self.mesh, "DG", 0, (self.I,self.I+self.K))
          self.dGammaidus_i = fd.Function(V_dGammaidus_i)
          self.dGammaidus_i.rename("PhaseSourceDerivativesFunction", "dGammaidus_i")
          self.dfunction_map[self.dGammaidus_i] = self.e_dGammaidus_i

          V_Gammaik_i = fd.VectorFunctionSpace(self.mesh, "DG", 0, self.K)
          self.Gammaik_i = fd.Function(V_Gammaik_i)
          self.Gammaik_i.rename("ComponentSourcesFunction", "Gammaik_i")
          self.function_map[self.Gammaik_i] = self.e_Gammaik_i

          V_dGammaikdus_i = fd.TensorFunctionSpace(self.mesh, "DG", 0, (self.K,self.I+self.K))
          self.dGammaikdus_i = fd.Function(V_dGammaikdus_i)
          self.dGammaikdus_i.rename("ComponentSourceDerivativesFunction", "dGammaikdus_i")
          self.dfunction_map[self.dGammaikdus_i] = self.e_dGammaikdus_i
        else:
          self.rhoi_i        = self.e_rhoi_i
          self.drhoidus_i    = self.e_drhoidus_i
          self.Gammai_i      = self.e_Gammai_i
          self.dGammaidus_i  = self.e_dGammaidus_i
          self.Gammik_i      = self.e_Gammaik_i
          self.dGammaikdus_i = self.e_dGammaikdus_i
        
        v_i = sum([self.mi_i[i]/self.rhoi_i[i] for i in range(self.I)])
        self.F = 0
        self.G = 0
        sKi = 0
        for i in range(self.I):
            self.F += self.mi_t[i]*self.mi_dot[i]*fd.dx
            self.G += self.mi_t[i]*self.Da_i*self.rho0_i*self.Gammai_i[i]*v_i*fd.dx
            Ki = self.cik_i[i].ufl_shape[0]
            for k in range(Ki):
                GikcGi_i = self.Gammaik_i[sKi+k] - self.cik_i[i][k]*self.Gammai_i[i]
                self.F += self.cik_t[i][k]*self.cik_dot[i][k]*fd.dx
                self.G += self.cik_t[i][k]*self.Da_i*self.rho0_i*GikcGi_i*v_i/(self.mi_i[i] + self.eps_i)*fd.dx
            sKi += Ki

        cd = {
              self.Gammai_i : self.dGammaidus_i, \
              self.Gammaik_i : self.dGammaikdus_i, \
              self.rhoi_i : self.drhoidus_i, \
             }
        self.JF = self.a_i*fd.derivative(self.F, self.us_dot, self.us_a, cd) + fd.derivative(self.F, self.us_i, self.us_a, cd)
        self.JG = fd.derivative(self.G, self.us_i, self.us_a, cd)

    def setup_mats(self):
        from firedrake.assemble import allocate_matrix, OneFormAssembler, TwoFormAssembler
        self.JFmat = allocate_matrix(self.JF, mat_type='aij')
        self.Fres = fd.Function(self.V)
        self.JGmat = allocate_matrix(self.JG, mat_type='aij')
        self.Gres = fd.Function(self.V)

        with self.Fres.dat.vec_wo as f:
          self.ts.setIFunction(self.FormIFunction, f)
        self.ts.setIJacobian(self.FormIJacobian, self.JFmat.petscmat, self.JFmat.petscmat)
        with self.Gres.dat.vec_wo as g:
          self.ts.setRHSFunction(self.FormRHSFunction, g)
        self.ts.setRHSJacobian(self.FormRHSJacobian, self.JGmat.petscmat, self.JGmat.petscmat)

        self.assemble_JF = TwoFormAssembler(self.JF, self.JFmat).assemble
        self.assemble_F  = OneFormAssembler(self.F,  self.Fres ).assemble
        self.assemble_JG = TwoFormAssembler(self.JG, self.JGmat).assemble
        self.assemble_G  = OneFormAssembler(self.G,  self.Gres ).assemble

    def setup_bounds(self):
        ub = self.us_i.dof_dset.layout_vec.duplicate()
        ub.setArray(np.ones(self.I+self.K))
        lb = self.us_i.dof_dset.layout_vec.duplicate()
        lb.setArray(np.zeros(self.I+self.K))
        self.snes.setVariableBounds(lb, ub)

    def set_initial_params(self,T,p,mi0,Cik0,**kwargs):
        super().set_initial_params(T,p,mi0,Cik0,**kwargs)
        self.T_i.assign(self.T)
        self.p_i.assign(self.p)
        self.Da_i.assign(self.Da)
        self.eps_i.assign(self.eps)
        self.rho0_i.assign(self.rho0)
        self.us_i.vector().set_local(np.concatenate((self.mi0,self.Cik0)))

    def update_functions(self):
        for f,e in self.function_map.items():
            f.vector().set_local(e.evaluate(np.ones(1)*0.5))

    def update_dfunctions(self):
        for f,e in self.dfunction_map.items():
            f.vector().set_local(e.evaluate(np.ones(1)*0.5))

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

        self.us_i.vector().set_local(np.concatenate((mi0, Cik0)))
        self.ts.setType(method)
        self.ts.setMaxTime(end)
        self.ts.setTolerances(rtol=rtol,atol=atol)
        self.ts.setMaxSteps(maxsteps)
        self.snes.setTolerances(rtol=snes_rtol, atol=snes_atol, max_it=snes_maxit)
        x = self.us_i.dof_dset.layout_vec.duplicate()
        #try:
        so = io.StringIO()
        se = io.StringIO()
        #with redirect_stdout(so), redirect_stderr(se):
        with self.us_i.dat.vec as v:
          v.copy(x)
          tic = time.perf_counter()
          self.ts.solve(x)
          toc = time.perf_counter()
          self.stime = toc-tic
          x.copy(v)
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
        with self.us_i.dat.vec_wo as v:
          u.copy(v)
        #y = self.us_i.at(0.5*np.ones(1))
        y = self.us_i.vector().array()
        self.sol.update(i,t,y)
        
    def FormIFunction(self, ts, t, u, udot, f):
        print('  In FormIFunction')
        
        with self.us_i.dat.vec_wo as v:
          u.copy(v)

        import ipdb; ipdb.set_trace()
        self.assemble_F()

        with self.Fres.dat.vec_ro as v:
          v.copy(f)

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

        with self.us_dot.dat.vec_wo as v:
          udot.copy(v)
        with self.us_i.dat.vec_wo as v:
          u.copy(v)
        
        self.assemble_JF()
        
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

        with self.us_i.dat.vec_wo as v:
          u.copy(v)
        
        self.update_functions()
        
        import ipdb; ipdb.set_trace()
        self.assemble_G()

        with self.Gres.dat.vec_ro as v:
          v.copy(g)

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

        with self.us_i.dat.vec_wo as v:
          u.copy(v)
        
        self.update_functions()
        self.update_dfunctions()

        self.assemble_JG()
        
        if self.print_norms:
            print('    FormRHSJacobian: 2-norm u = {}'.format(u.norm(norm_type=1)))
            print('    FormRHSJacobian: inf-norm u = {}'.format(u.norm(norm_type=3)))
            print('    FormRHSJacobian: Frobenius norm A = {}'.format(gmat.norm(norm_type=2)))
            print('    FormRHSJacobian: inf-norm A = {}'.format(gmat.norm(norm_type=3)))
            print('    FormRHSJacobian: Frobenius norm B = {}'.format(pgmat.norm(norm_type=2)))
            print('    FormRHSJacobian: inf-norm B = {}'.format(pgmat.norm(norm_type=3)))
        return True

if not FIREDRAKE_AVAIL:
    class FiredrakePDReactiveODE(PETScPDReactiveODE):
        def __init__(self,rxn):
            raise Exception("Firedrake not available")

