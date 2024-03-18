import io
import numpy as np
from scipy.integrate import solve_ivp
from contextlib import redirect_stdout, redirect_stderr
import time
import inspect

#from ..base import *

import dolfinx
import basix
import ufl
from petsc4py import PETSc
from mpi4py import MPI

class TCGReaction:
  rxn = None
  I = None
  Kis = None
  K = None
  J = None
  mi = None
  cik = None

  def __init__(self, rxn):
    self.rxn = rxn
    self.I = len(self.rxn.phases())
    self.Kis = [len(self.rxn.phases()[i].endmembers()) for i in range(self.I)]
    self.K = sum(self.Kis)
    self.J = len(self.rxn.nu())
    self.mi = np.zeros(self.I)
    self.cik = self.rxn.zero_C()

class TCGCoefficients(TCGReaction):
  mi_funcs = None
  cik_funcs = None
  T_func = None
  p_func = None
  eps_const = None

  rhoi_i = None
  Gammai_i = None
  Gammaik_i = None
  drhoidus_i = None
  dGammaidus_i = None
  dGammaikdus_i = None

  eps = None
  T = None 
  p = None

  T_values = None
  p_values = None
  mi_values = None
  cik_values = None

  mesh = None
  element = None #?

  default_x = None
  default_cells = None

  params_evald = False

  def __init__(self, mi_funcs, cik_funcs, T_func, p_func, eps_const, mesh, element, tcgrxn):
    super().__init__(tcgrxn.rxn)
    
    self.mi_funcs = mi_funcs
    self.cik_funcs = cik_funcs
    self.T_func = T_func
    self.p_func = p_func
    self.eps_const = eps_const
    self.mesh = mesh
    self.element = element

    rhoi_e = ufl.VectorElement(self.element, dim=self.I)
    rhoi_V = dolfinx.fem.FunctionSpace(self.mesh, rhoi_e)
    self.rhoi_i = dolfinx.fem.Function(rhoi_V)

    drhoidus_e = ufl.TensorElement(self.element, shape=(self.I, self.I+self.K))
    drhoidus_V = dolfinx.fem.FunctionSpace(self.mesh, drhoidus_e)
    self.drhoidus_i = dolfinx.fem.Function(drhoidus_V)

    Gammai_e = ufl.VectorElement(self.element, dim=self.I)
    Gammai_V = dolfinx.fem.FunctionSpace(self.mesh, Gammai_e)
    self.Gammai_i = dolfinx.fem.Function(Gammai_V)

    dGammaidus_e = ufl.TensorElement(self.element, shape=(self.I, self.I+self.K))
    dGammaidus_V = dolfinx.fem.FunctionSpace(self.mesh, dGammaidus_e)
    self.dGammaidus_i = dolfinx.fem.Function(dGammaidus_V)

    Gammaik_e = ufl.VectorElement(self.element, dim=self.K)
    Gammaik_V = dolfinx.fem.FunctionSpace(self.mesh, Gammaik_e)
    self.Gammaik_i = dolfinx.fem.Function(Gammaik_V)

    dGammaikdus_e = ufl.TensorElement(self.element, shape=(self.K, self.I+self.K))
    dGammaikdus_V = dolfinx.fem.FunctionSpace(self.mesh, dGammaikdus_e)
    self.dGammaikdus_i = dolfinx.fem.Function(dGammaikdus_V)

    tdim = self.mesh.topology.dim
    num_cells = self.mesh.topology.index_map(tdim).size_local + self.mesh.topology.index_map(tdim).num_ghosts

    gdim = self.mesh.geometry.dim
    x_dofmap = self.mesh.geometry.dofmap
    num_dofs_g = len(x_dofmap.links(0))
    x_g = self.mesh.geometry.x

    cmap = self.mesh.geometry.cmap
    # need to convert to a basix element to get the tabulate function
    # FIXME: assumed Lagange and continuous and no attempt to include variants!
    basix_cmap = basix.create_element(basix.finite_element.string_to_family("Lagrange", self.mesh.topology.cell_name()), 
                                      basix.cell.string_to_type(self.mesh.topology.cell_name()), cmap.degree, False)

    # need to access the element through a functionspace here as it changes type and allows access to interpolation points
    X = rhoi_V.element.interpolation_points
    phi = basix_cmap.tabulate(0, X)[0,:,:,0]

    coordinate_dofs = np.zeros((num_dofs_g, gdim))

    self.default_x = np.zeros((3, num_cells * X.shape[0]))
    self.default_cells = np.empty(num_cells * X.shape[0], dtype=int)

    for c in range(num_cells):
      x_dofs = x_dofmap.links(c)
      for i in range(len(x_dofs)): coordinate_dofs[i][:] = x_g[x_dofs[i]][:gdim]
      for p in range(X.shape[0]):
        for j in range(gdim):
          self.default_x[j, c*X.shape[0]+p] = sum([phi[p,k] * coordinate_dofs[k,j] for k in range(num_dofs_g)])
        self.default_cells[c * X.shape[0] + p] = c

    self.T_values = np.empty((self.default_x.shape[1], 1), dtype=PETSc.ScalarType)
    self.p_values = np.empty((self.default_x.shape[1], 1), dtype=PETSc.ScalarType)
    self.mi_values = [np.empty((self.default_x.shape[1], 1), dtype=PETSc.ScalarType) for i in range(self.I)]
    self.cik_values = [np.empty((self.default_x.shape[1], self.Kis[i]), dtype=PETSc.ScalarType) for i in range(self.I)]

  def interpolate_all_coeffs(self):
    self.interpolate_coeffs()
    self.params_evald = True
    self.interpolate_derivative_coeffs()

  def interpolate_coeffs(self):
    self.eval_params(self.default_x)
    self.params_evald = True
    self.rhoi_i.interpolate(self.rhoi)
    self.Gammai_i.interpolate(self.Gammai)
    self.Gammaik_i.interpolate(self.Gammaik)
    self.params_evald = False

  def interpolate_derivative_coeffs(self):
    self.eval_params(self.default_x)
    self.params_evald = True
    self.drhoidus_i.interpolate(self.drhoidus)
    self.dGammaidus_i.interpolate(self.dGammaidus)
    self.dGammaikdus_i.interpolate(self.dGammaikdus)
    self.params_evald = False

  def get_interpolation_cells(self, x):
    if x.shape == self.default_x.shape:
      return self.default_cells
    else:
      tdim = self.mesh.topology.dim
      num_cells = self.mesh.topology.index_map(tdim).size_local + self.mesh.topology.index_map(tdim).num_ghosts
      cell_indices = np.arange(num_cells, dtype=np.int32)
      midpoint_tree = dolfinx.geometry.create_midpoint_tree(self.mesh, tdim, cell_indices)
      tree = dolfinx.geometry.BoundingBoxTree(self.mesh, tdim)
      return dolfinx.geometry.compute_closest_entity(tree, midpoint_tree, self.mesh, x.T)

  def eval_params(self, x):
    if not self.params_evald:
      cells = self.get_interpolation_cells(x)
      xT = x.transpose()

      self.eval_eps()
      self.eval_T(xT, cells)
      self.eval_p(xT, cells)
      self.eval_mi(xT, cells)
      self.eval_cik(xT, cells)

  def eval_eps(self):
    self.eps = PETSc.ScalarType(self.eps_const.value)

  def eval_T(self, xT, cells):
    if self.T_values.shape != (xT.shape[0], 1):
      self.T_values = np.empty((xT.shape[0], 1), dtype=PETSc.ScalarType)
    self.T_values = self.T_func.eval(xT, cells, self.T_values).reshape(xT.shape[0], 1)
    
  def eval_p(self, xT, cells):
    if self.p_values.shape != (xT.shape[0], 1):
      self.p_values = np.empty((xT.shape[0], 1), dtype=PETSc.ScalarType)
    self.p_values = self.p_func.eval(xT, cells, self.p_values).reshape(xT.shape[0], 1)
    
  def eval_mi(self, xT, cells):
    for i in range(self.I):
      if self.mi_values[i].shape != (xT.shape[0], 1):
        self.mi_values[i] = np.empty((xT.shape[0], 1), dtype=PETSc.ScalarType)
      self.mi_values[i] = self.mi_funcs[i].eval(xT, cells, self.mi_values[i]).reshape(xT.shape[0], 1)
    
  def eval_cik(self, xT, cells):
    for i in range(self.I):
      if self.cik_values[i].shape != (xT.shape[0], self.Kis[i]):
        self.cik_values[i] = np.empty((xT.shape[0], self.Kis[i]), dtype=PETSc.ScalarType)
      self.cik_values[i] = self.cik_funcs[i].eval(xT, cells, self.cik_values[i]).reshape(xT.shape[0], self.Kis[i])

  def update_params_x(self, xi):
    self.update_T_x(xi)
    self.update_p_x(xi)
    self.update_mi_x(xi)
    self.update_cik_x(xi)
    
  def update_T_x(self, xi):
    self.T = self.T_values[xi,0]

  def update_p_x(self, xi):
    self.p = self.p_values[xi,0]

  def update_mi_x(self, xi):
    self.mi = np.asarray([self.mi_values[i][xi,0] for i in range(self.I)])

  def update_cik_x(self, xi):
    for i in range(self.I):
      self.cik[i] = np.clip(self.cik_values[i][xi,:], self.eps, 1.-self.eps)
      self.cik[i] /= np.sum(self.cik[i])

  def rho(self, x):
    x = np.asarray(x)
    rhoi = self.rhoi(x) # eval_params taken care of here if necessary
    xvals = np.empty(x.shape[1])
    for xi in range(x.shape[1]):
      self.update_params_x(xi)
      xvals[xi] = 1./sum(self.mi/rhoi[:,xi])
    return xvals

  def rhoi(self, x):
    self.eval_params(x)
    xvals = np.empty((self.I, x.shape[1]))
    for xi in range(x.shape[1]):
      self.update_params_x(xi)
      xvals[:,xi] = np.asarray(self.rxn.rho(self.T, self.p, self.cik))
    return xvals
        
  def drhoidus(self, x):
    self.eval_params(x)
    xvals = np.zeros((self.I*(self.I+self.K), x.shape[1]))
    for xi in range(x.shape[1]):
      self.update_params_x(xi)
      for i in range(self.I):
        xvals[i*(self.I+self.K) + self.I + sum(self.Kis[:i]):i*(self.I+self.K) + self.I + sum(self.Kis[:i+1]), xi] = \
                                  self.rxn.phases()[i].drho_dc(self.T, self.p, self.cik[i])
    return xvals
  
  def Gammai(self, x):
    self.eval_params(x)
    xvals = np.empty((self.I, x.shape[1]))
    for xi in range(x.shape[1]):
      self.update_params_x(xi)
      xvals[:,xi] = np.asarray(self.rxn.Gamma_i(self.T,self.p,self.cik,self.mi))
    return xvals
  
  def dGammaidus(self, x):
    self.eval_params(x)
    xvals = np.empty((self.I*(self.I+self.K), x.shape[1]))
    for xi in range(x.shape[1]):
      self.update_params_x(xi)
      
      dGammaidmi = np.asarray(self.rxn.dGamma_i_dPhi(self.T,self.p,self.cik,self.mi))
      dGammaidcik = self.rxn.dGamma_i_dC(self.T,self.p,self.cik,self.mi)

      for i in range(self.I):
        xvals[i*(self.I+self.K):i*(self.I+self.K)+self.I,xi] = dGammaidmi[i][:]
        for di in range(self.I):
          xvals[i*(self.I+self.K) + self.I + sum(self.Kis[:di]):i*(self.I+self.K) + self.I + sum(self.Kis[:di+1]), xi] = \
                   dGammaidcik[i][di][:]
    return xvals
  
  def Gammaik(self, x):
    self.eval_params(x)
    xvals = np.empty((self.K, x.shape[1]))
    for xi in range(x.shape[1]):
      self.update_params_x(xi)

      Gammaik = self.rxn.Gamma_ik(self.T, self.p, self.cik, self.mi)

      for i in range(self.I):
        xvals[sum(self.Kis[:i]):sum(self.Kis[:i+1]),xi] = Gammaik[i][:]
    return xvals
  
  def dGammaikdus(self, x):
    self.eval_params(x)
    xvals = np.empty((self.K*(self.I+self.K), x.shape[1]))
    for xi in range(x.shape[1]):
      self.update_params_x(xi)
      
      for i in range(self.I):
        dGammaikdmii = self.rxn.dGamma_ik_dPhi(self.T,self.p,self.cik,self.mi,i)
        for k in range(self.Kis[i]):
          xvals[(sum(self.Kis[:i]) + k)*(self.I + self.K):(sum(self.Kis[:i]) + k)*(self.I + self.K) + self.I, xi] = \
                dGammaikdmii[k][:]
          for di in range(self.I):
            for dk in range(self.Kis[di]):
              xvals[(sum(self.Kis[:i]) + k)*(self.I + self.K) + self.I + sum(self.Kis[:di]) + dk, xi] = \
                  self.rxn.dGamma_ik_dC(self.T,self.p,self.cik,self.mi,i,k,di,dk)
    return xvals

class PDReactiveODE:
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

  sol = None

  # default values
  i0 = 0
  ci0k0 = 0.9
  T = 1673.
  p = 150000.
  eps = 1.e-2
  Da = 1.0
  x0 = [[0.5], [0.5], [0.5]]

  def __init__(self,rxn,nx=1):
    self.comm = MPI.COMM_WORLD

    # reaction object
    self.tcgrxn = TCGReaction(rxn)

    # mesh
    self.mesh = dolfinx.mesh.create_interval(self.comm, nx, [0.,1.], dolfinx.mesh.GhostMode.none)

    # elements and system function space
    mi_es = [ufl.FiniteElement("DG", self.mesh.ufl_cell(), 0) for i in range(self.tcgrxn.I)]
    cik_es = [ufl.VectorElement(ufl.FiniteElement("DG", self.mesh.ufl_cell(), 0), dim=self.tcgrxn.Kis[i]) for i in range(self.tcgrxn.I)]
    us_e = ufl.MixedElement(mi_es + cik_es)
    self.us_V = dolfinx.fem.FunctionSpace(self.mesh, us_e)

    # solution vector
    self.us_i = dolfinx.fem.Function(self.us_V)
    self.us_i.name = "us_i"
    us_is = ufl.split(self.us_i)
    mi_is = us_is[:self.tcgrxn.I]
    cik_is = us_is[self.tcgrxn.I:]

    self.mi_i  = [self.us_i.sub(i) for i in range(self.tcgrxn.I)]
    for i in range(self.tcgrxn.I): self.mi_i[i].name = self.tcgrxn.rxn.phases()[i].name()+"MassFraction"
    self.cik_i = [self.us_i.sub(i) for i in range(self.tcgrxn.I,2*self.tcgrxn.I)]
    for i in range(self.tcgrxn.I): self.cik_i[i].name = self.tcgrxn.rxn.phases()[i].name()+"ComponentMassFractions"

    # temperature
    T_e = ufl.FiniteElement("DG", self.mesh.ufl_cell(), 0)
    T_V = dolfinx.fem.FunctionSpace(self.mesh, T_e)
    self.T_i = dolfinx.fem.Function(T_V)
    self.T_i.name = "Temperature"

    # pressure
    p_e = ufl.FiniteElement("DG", self.mesh.ufl_cell(), 0)
    p_V = dolfinx.fem.FunctionSpace(self.mesh, p_e)
    self.p_i = dolfinx.fem.Function(p_V)
    self.p_i.name = "Pressure"

    # eps is necessary to initialize the coeffs
    self.eps_i  = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(self.eps))

    # set up the coeffs class
    coeffs_element = ufl.FiniteElement("DG", self.mesh.ufl_cell(), 0)
    self.coeffs = TCGCoefficients(self.mi_i, self.cik_i, 
                                  self.T_i, self.p_i, 
                                  self.eps_i, self.mesh,
                                  coeffs_element, self.tcgrxn)

    # the rest of the constants
    self.Da_i   = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(self.Da))
    self.a_i    = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(666.e6))
    self.rho0_i = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(666.e6)) # initial condition needs to be set first

    # time-derivatives
    self.us_dot = dolfinx.fem.Function(self.us_V)
    self.us_dot.name = "us_dot"
    self.us_dot.x.set(0.0)
    us_dots = ufl.split(self.us_dot)
    mi_dots = us_dots[:self.tcgrxn.I]
    cik_dots = us_dots[self.tcgrxn.I:]

    # function views into time-derivatives
    self.mi_dot  = [self.us_dot.sub(i) for i in range(self.tcgrxn.I)]
    for i in range(self.tcgrxn.I): self.mi_dot[i].name = self.tcgrxn.rxn.phases()[i].name()+"MTimeDerivative"
    self.cik_dot = [self.us_dot.sub(i) for i in range(self.tcgrxn.I,2*self.tcgrxn.I)]
    for i in range(self.tcgrxn.I): self.cik_dot[i].name = self.tcgrxn.rxn.phases()[i].name()+"CTimeDerivative"

    # test functions
    self.us_t  = ufl.TestFunction(self.us_V)
    us_ts = ufl.split(self.us_t)
    mi_ts  = us_ts[:self.tcgrxn.I]
    cik_ts = us_ts[self.tcgrxn.I:]

    # trial function
    us_a  = ufl.TrialFunction(self.us_V)

    # forms!
    v_i = sum([mi_is[i]/self.coeffs.rhoi_i[i] for i in range(self.tcgrxn.I)])
    F = 0
    G = 0
    sKi = 0
    for i in range(self.tcgrxn.I):
      F += mi_ts[i]*mi_dots[i]*ufl.dx
      G += mi_ts[i]*self.Da_i*self.rho0_i*self.coeffs.Gammai_i[i]*v_i*ufl.dx
      Ki = cik_is[i].ufl_shape[0]
      for k in range(Ki):
        GikcGi_i = self.coeffs.Gammaik_i[sKi+k] - cik_is[i][k]*self.coeffs.Gammai_i[i]
        F += cik_ts[i][k]*cik_dots[i][k]*ufl.dx
        G += cik_ts[i][k]*self.Da_i*self.rho0_i*GikcGi_i*v_i/(mi_is[i] + self.eps_i)*ufl.dx
      sKi += Ki

    # coefficient derivatives
    cd = {
          self.coeffs.rhoi_i    : self.coeffs.drhoidus_i,    \
          self.coeffs.Gammai_i  : self.coeffs.dGammaidus_i,  \
          self.coeffs.Gammaik_i : self.coeffs.dGammaikdus_i, \
         }
    JF = self.a_i*ufl.derivative(F, self.us_dot, us_a, cd) + ufl.derivative(F, self.us_i, us_a, cd)
    JG = ufl.derivative(G, self.us_i, us_a, cd)

    self.F  = dolfinx.fem.form(F)
    self.G  = dolfinx.fem.form(G)
    self.JF = dolfinx.fem.form(JF)
    self.JG = dolfinx.fem.form(JG)

  def initialize(self):
    # apply initial conditions
    for i in range(self.tcgrxn.I): 
      self.mi_i[i].interpolate(lambda x: np.ones((1,x.shape[1])) if i==self.i0 else np.zeros((1,x.shape[1])))
    def ci0k(x):
        ci0k = np.empty((self.tcgrxn.Kis[i],1))
        if self.tcgrxn.Kis[i]==1:
            ci0k[0] = 1
        else:
            ci0k[0] = self.ci0k0
        for k in range(1,self.tcgrxn.Kis[i]): ci0k[k] = (1.-self.ci0k0)/(self.tcgrxn.Kis[i]-1.)
        return np.tile(ci0k, (1,x.shape[1]))
    for i in range(self.tcgrxn.I): self.cik_i[i].interpolate(ci0k)

    # re-initialize constants but most important is rho0, which cannot be set until initial conditions are known
    self.T_i.x.set(self.T)
    self.p_i.x.set(self.p)
    self.rho0_i.value = PETSc.ScalarType(self.coeffs.rho(self.x0)[0])
    self.eps_i.value  = PETSc.ScalarType(self.eps)
    self.Da_i.value   = PETSc.ScalarType(self.Da)

    self.sol = self.Solution(self.tcgrxn.I+self.tcgrxn.K)

  def work_vector(self):
    return dolfinx.cpp.la.petsc.create_vector_wrap(self.us_i.x)

  def bounds(self):
    dofmap = self.us_V.dofmap
    lb = dolfinx.la.create_petsc_vector(dofmap.index_map, dofmap.index_map_bs)
    lb.zeroEntries()
    ub = dolfinx.la.create_petsc_vector(dofmap.index_map, dofmap.index_map_bs)
    ub.set(1.0)
    return lb, ub

  def interpolate_ifunction_coeffs(self):
    pass

  def interpolate_rhsfunction_coeffs(self):
    self.coeffs.interpolate_coeffs()

  def interpolate_ijacobian_coeffs(self):
    pass

  def interpolate_rhsjacobian_coeffs(self):
    self.coeffs.interpolate_all_coeffs()

class TSProblem:

  # default parameters
  print_norms = False

  def __init__(self, problem):

    self.problem = problem
    self.comm = self.problem.comm

    self.ts = PETSc.TS().create(comm=self.comm)
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
    self.pc.setFactorSolverType('umfpack')

    self.ts.setFromOptions()
    
    self.Fres = dolfinx.fem.petsc.create_vector(self.problem.F)
    self.ts.setIFunction(self.FormIFunction, self.Fres)

    self.JFmat = dolfinx.fem.petsc.create_matrix(self.problem.JF)
    self.ts.setIJacobian(self.FormIJacobian, self.JFmat, self.JFmat)

    self.Gres = dolfinx.fem.petsc.create_vector(self.problem.G)
    self.ts.setRHSFunction(self.FormRHSFunction, self.Gres)

    self.JGmat = dolfinx.fem.petsc.create_matrix(self.problem.JG)
    self.ts.setRHSJacobian(self.FormRHSJacobian, self.JGmat, self.JGmat)

    lb, ub = self.problem.bounds()
    self.snes.setVariableBounds(lb, ub)

  def solve(self,**kwargs):
    method = kwargs.get('ts_type', 'bdf')
    rtol   = kwargs.get('ts_rtol', 1.e-5)
    atol   = kwargs.get('ts_atol', 1.e-9)
    end    = kwargs.get('ts_max_time', 1.0)
    snes_rtol = kwargs.get('snes_rtol', 1.e-6)
    snes_atol = kwargs.get('snes_atol', 1.e-6)
    snes_maxit = kwargs.get('snes_maxit', 10)
    maxsteps = kwargs.get('maxsteps', 100000000)

    self.print_norms = kwargs.get('print_norms', False)

    self.ts.setType(method)
    self.ts.setMaxTime(end)
    self.ts.setTolerances(rtol=rtol,atol=atol)
    self.ts.setMaxSteps(maxsteps)
    self.snes.setTolerances(rtol=snes_rtol, atol=snes_atol, max_it=snes_maxit)

    self.problem.initialize()
    x = self.problem.work_vector()

    try:
      so = io.StringIO()
      se = io.StringIO()
      #with redirect_stdout(so), redirect_stderr(se):
      tic = time.perf_counter()
      self.ts.solve(x)
      reason = self.ts.getConvergedReason()
      print("  TSConvergedReason: ", reason)
      toc = time.perf_counter()
      self.stime = toc-tic
      #  self.stdout = so.getvalue()
      #  self.stderr = se.getvalue()
      flag = self.problem.tcgrxn.rxn.check_coder_error()
      if flag!=0:
        self.sol = None
        self.excstr = repr(flag)
        self.rxn.reset_coder_error()
    except Exception as e:
      self.sol = None
      self.excstr = repr(e)

    dolfinx.common.list_timings(MPI.COMM_WORLD, [dolfinx.common.TimingType.wall])

#  def plot(self):
#    fig = plt.figure(figsize=(12,24))
#    axes = [fig.add_subplot(3,2,i+1) for i in  range(6)]
#    
#    self.plot_axes(axes)
#    self.apply_labels(axes)
#
#    fig.suptitle('$p$ = {:.1f} GPa, $T$ = {:.1f} K, $\\rho_0$ = {:.1f} kg/m$^3$, eps = {:.2e}'.format(\
#                 Bar2GPa(self.problem.p),self.problem.T, rho2kgpm3(float(self.problem.rho0_i.value)),self.eps),y=0.9)
#    
#    plt.show()
#  
#  def apply_labels(self,axes):
#      
#    plabels = [self.problem.tcgrxn.rxn.phases()[i].abbrev() for i in range(self.problem.tcgrxn.I)]
#    elabels = [self.problem.tcgrxn.rxn.phases()[i].endmembers()[k].formula()+'_('+self.problem.tcgrxn.rxn.phases()[i].abbrev()+')' \
#                                       for i in range(self.problem.tcgrxn.I) for k in range(self.problem.tcgrxn.Kis[i])]
#
#    axes[0].set_ylabel('$\\phi_i$')
#    labels = [line.get_label()  for line in  axes[0].get_lines()[len(plabels):]]
#    axes[0].legend(plabels+labels,loc='best')
##         axes[0].set_ylim([-0.05,1.05])
#    
#    axes[1].set_ylabel('$\Gamma_i$')
#    
#    axes[2].set_ylabel('$m_i$')
#    
#    axes[3].set_ylabel('$\\rho$ kg/m$^3$')
#    axes[3].legend(loc='best')
#    
#    axes[4].set_ylabel('$C_i^k$')
#    axes[4].set_ylim([-0.05,1.05])
#    axes[4].legend(elabels, loc='best')
#    
#    axes[5].set_ylabel('$\Gamma_i^k$')
#    
#    for axis in axes:
#      axis.grid()
#      axis.set_xlabel('$t^*$')
#              
#  def plot_axes(self, axes, ls='-', label=None):
#    assert(self.sol is not None)
#    assert(len(axes)==6)
#    
#    sol = self.sol
#    
#    t = self.sol.t
#    
#    T = self.T
#    p = self.p
#    
#    mi  = self.sol.y[:self.problem.tcgrxn.I].T
#    cik = self.sol.y[self.problem.tcgrxn.I:self.problem.tcgrxn.I+self.problem.tcgrxn.K].T
#
#    indices_eval = list(range(len(t)))
#    gamma_i  = np.zeros(mi.shape)
#    gamma_ik = np.zeros(cik.shape)
#    rhoi     = gamma_i.copy()
#    rho      = np.zeros(t.shape)
#    c_ik     = np.zeros(cik.shape)
#
#    for i in range(len(t)):
#      C = self.reshapeC(Cik[i])
#      Cs = self.regularizeC(C)
#      gamma_i[i,:]  = self.rxn.Gamma_i(T,p,Cs,mi[i])
#      gamma_ik[i,:] = [gik for gi in self.rxn.Gamma_ik(T,p,Cs,mi[i]) for gik in gi]
#      rhoi[i]       = self.rxn.rho(T, p, Cs)
#      v             = mi[i]/rhoi[i]
#      rho[i]        = 1./v.sum()
#      C_ik[i,:]     = [cik for ci in Cs for cik in ci]
#
#    #for i in range(self.I): axes[0].plot(t,rho*mi[:,i]/rhoi[:,i],color=colors[i],ls=ls)
#    for i in range(self.I): axes[0].plot(t,rho*mi[:,i]/rhoi[:,i],ls=ls)
#    mlabel='$\\sum_i\\phi_i$'
#    if label is not None: mlabel=mlabel+' '+label
#    sumphiline = axes[0].plot(t,rho*(mi/rhoi).sum(axis=-1),color='k',ls='-',label=mlabel)
#    
#    #for i in range(self.I): axes[1].plot(t,gamma_i[:,i],color=colors[i],ls=ls)
#    for i in range(self.I): axes[1].plot(t,gamma_i[:,i],ls=ls)
#
##         for i in range(self.I): 
##             for k in range(self.Kis[i]):
##                 axes[4].plot(t,mi[:,i]*Cik[:,sum(self.Kis[:i])+k],ls=ls,color=colors[sum(self.Kis[:i])+k])
#
#    #for i in range(self.I): axes[2].plot(t,mi[:,i],color=colors[i],ls=ls)
#    for i in range(self.I): axes[2].plot(t,mi[:,i],ls=ls)
#
#    #for i in range(self.I): axes[3].plot(t,rho2kgpm3(rhoi[:,i]),color=colors[i],ls=ls,label='_')
#    for i in range(self.I): axes[3].plot(t,rho2kgpm3(rhoi[:,i]),ls=ls,label='_')
#    mlabel='$\\rho$'
#    if label is not None: mlabel=mlabel+' '+label
#    rholine = axes[3].plot(t,rho2kgpm3(rho),'k',ls=ls,label=mlabel)
#
#    #for ik in range(self.K): axes[4].plot(t,Cik[:,ik],ls=ls,color=colors[len(colors)%(ik+1)])
#    for ik in range(self.K): axes[4].plot(t,Cik[:,ik],ls=ls)
#
#    #for ik in range(self.K): axes[5].plot(t,gamma_ik[:,ik],ls=ls,color=colors[len(colors)%(ik+1)])
#    for ik in range(self.K): axes[5].plot(t,gamma_ik[:,ik],ls=ls)

  def FormIFunction(self, ts, t, u, udot, f):
    print('  In FormIFunction')

    uu    = u.getArray(readonly=True)
    uudot = udot.getArray(readonly=True)

    self.problem.us_i.x.array[:]   = uu[:]
    self.problem.us_dot.x.array[:] = uudot[:]

    self.problem.interpolate_ifunction_coeffs()

    f.zeroEntries()
    dolfinx.fem.assemble_vector(f.getArray(), self.problem.F)
    
    if self.print_norms:
      print('    FormIFunction: 2-norm u = {}'.format(u.norm(norm_type=1)))
      print('    FormIFunction: inf-norm u = {}'.format(u.norm(norm_type=3)))
      print('    FormIFunction: 2-norm u_t = {}'.format(udot.norm(norm_type=1)))
      print('    FormIFunction: inf-norm u_t = {}'.format(udot.norm(norm_type=3)))
      print('    FormIFunction: 2-norm f = {}'.format(f.norm(norm_type=1)))
      print('    FormIFunction: inf-norm f = {}'.format(f.norm(norm_type=3)))
      
  def FormRHSFunction(self, ts, t, u, g):
    print('  In FormRHSFunction')

    uu = u.getArray(readonly=True)
    self.problem.us_i.x.array[:] = uu[:]

    self.problem.interpolate_rhsfunction_coeffs()

    g.zeroEntries()
    dolfinx.fem.assemble_vector(g.getArray(), self.problem.G)

    if self.print_norms:
      print('    FormRHSFunction: 2-norm u = {}'.format(u.norm(norm_type=1)))
      print('    FormRHSFunction: inf-norm u = {}'.format(u.norm(norm_type=3)))
      print('    FormRHSFunction: 2-norm g = {}'.format(g.norm(norm_type=1)))
      print('    FormRHSFunction: inf-norm g = {}'.format(g.norm(norm_type=3)))
      
  def FormIJacobian(self, ts, t, u, udot, shift, fmat, pfmat):
    print('  In FormIJacobian')
    print("    a (shift) = {}".format(shift))

    self.problem.a_i.value = shift

    uu    = u.getArray(readonly=True)
    uudot = udot.getArray(readonly=True)

    self.problem.us_i.x.array[:]   = uu[:]
    self.problem.us_dot.x.array[:] = uudot[:]

    self.problem.interpolate_ijacobian_coeffs()

    fmat.zeroEntries()
    dolfinx.fem.petsc.assemble_matrix(fmat, self.problem.JF)
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

  def FormRHSJacobian(self, ts, t, u, gmat, pgmat):
    print('  In FormRHSJacobian')

    uu    = u.getArray(readonly=True)

    self.problem.us_i.x.array[:]   = uu[:]

    self.problem.interpolate_rhsjacobian_coeffs()

    gmat.zeroEntries()
    dolfinx.fem.petsc.assemble_matrix(gmat, self.problem.JG)
    gmat.assemble()

    if self.print_norms:
      print('    FormRHSJacobian: 2-norm u = {}'.format(u.norm(norm_type=1)))
      print('    FormRHSJacobian: inf-norm u = {}'.format(u.norm(norm_type=3)))
      print('    FormRHSJacobian: Frobenius norm A = {}'.format(gmat.norm(norm_type=2)))
      print('    FormRHSJacobian: inf-norm A = {}'.format(gmat.norm(norm_type=3)))
      print('    FormRHSJacobian: Frobenius norm B = {}'.format(pgmat.norm(norm_type=2)))
      print('    FormRHSJacobian: inf-norm B = {}'.format(pgmat.norm(norm_type=3)))
    return True

  def SNESCustomMonitor(self,snes,n,norm):
    print('    {} SNES Function norm {}'.format(n,norm))
        
  def TSCustomMonitor(self,ts,i,t,u):
    dt = ts.getTimeStep()
    print('{} TS dt {} time {}'.format(i,dt,t))
    y = np.array(u[:].tolist())
    self.problem.sol.update(i,t,y)
        
