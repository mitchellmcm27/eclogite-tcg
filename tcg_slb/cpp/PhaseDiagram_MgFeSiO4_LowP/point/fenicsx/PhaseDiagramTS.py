import py_MgFeSiO4_all_slb_rx as tcgdb
rxn = tcgdb.MgFeSiO4_all_slb_rx()

from ufl import ( FiniteElement, VectorElement, TensorElement, MixedElement,
                  TestFunction, TrialFunction, 
                  Coefficient, Constant,
                  interval, split, dx, derivative )

I = len(rxn.phases())
Kis = [len(rxn.phases()[i].endmembers()) for i in range(I)]
K = sum(Kis)
J = len(rxn.nu())

# Element declaration
mi_es = [FiniteElement("DG", interval, 0) for i in range(I)]
cik_es = [VectorElement(FiniteElement("DG", interval, 0), dim=Kis[i]) for i in range(I)]
us_e = MixedElement(mi_es + cik_es)

# Test space declaration
us_t = TestFunction(us_e)
us_ts = split(us_t)
mi_t = us_ts[:I]
cik_t = us_ts[I:]

# Trial space declaration
us_a = TrialFunction(us_e)

# Last iteration value declaration for System: System
us_i = Coefficient(us_e)
us_is = split(us_i)
mi_i = us_is[:I]
cik_i = us_is[I:]

# Previous time-level value declaration for System: System
us_dot = Coefficient(us_e)
us_dots = split(us_dot)
mi_dot = us_dots[:I]
cik_dot = us_dots[I:]

# DamkoehlerNumber
Da = Constant(interval)
# a
a = Constant(interval)
# Epsilon
eps = Constant(interval)
# DensityScale
rho0 = Constant(interval)

coeff_e = FiniteElement("DG", interval, 0)

# PhaseDensities
rhoi_e = VectorElement(coeff_e, dim=I)
rhoi_i = Coefficient(rhoi_e)

# PhaseDensityDerivatives
drhoidus_e = TensorElement(coeff_e, shape=(I,I+K))
drhoidus_i = Coefficient(drhoidus_e)

# PhaseSources
Gammai_e = VectorElement(coeff_e, dim=I)
Gammai_i = Coefficient(Gammai_e)

# PhaseSourceDerivatives
dGammaidus_e = TensorElement(coeff_e, shape=(I,I+K))
dGammaidus_i = Coefficient(dGammaidus_e)

# ComponentSources
Gammaik_e = VectorElement(coeff_e, dim=K)
Gammaik_i = Coefficient(Gammaik_e)

# ComponentSourceDerivatives
dGammaikdus_e = TensorElement(coeff_e, shape=(K,I+K))
dGammaikdus_i = Coefficient(dGammaikdus_e)

v_i = sum([mi_i[i]/rhoi_i[i] for i in range(I)])

# Residual
F = 0
G = 0
sKi = 0
for i in range(I):
  F += mi_t[i]*mi_dot[i]*dx
  G += mi_t[i]*Da*rho0*Gammai_i[i]*v_i*dx
  Ki = cik_i[i].ufl_shape[0]
  for k in range(Ki):
    GikcGi_i = Gammaik_i[sKi+k] - cik_i[i][k]*Gammai_i[i]
    F += cik_t[i][k]*cik_dot[i][k]*dx
    G += cik_t[i][k]*Da*rho0*GikcGi_i*v_i/(mi_i[i] + eps)*dx
  sKi += Ki

# Jacobian
cd = { \
      Gammai_i  : dGammaidus_i,  \
      Gammaik_i : dGammaikdus_i, \
      rhoi_i    : drhoidus_i,    \
     }
JF = a*derivative(F, us_dot, us_a, cd) + derivative(F, us_i, us_a, cd)
JG = derivative(G, us_i, us_a, cd)

# Declare potentially non-default form names to be accessible
forms = [F, G, JF, JG]

