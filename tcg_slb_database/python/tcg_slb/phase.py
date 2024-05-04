import numpy
import sympy as sym
import pandas as pd
from thermocodegen.coder import coder

class SLBPhase:
    name   = None
    abbrev = None
    reference = None

    nc = None
    ns = None
    np = None
    
    endmember_names = None
    formula_str     = None
    conversion_strs = None
    test_strs       = None
    
    dsyms = None
    Wsarr = None
    m = None
    C = None
    
    # Landau parameters
    TC0 = None
    VD  = None
    SD  = None
    TC0syms = None
    VDsyms  = None
    SDsyms  = None
    
    Rsym = None
    
    G = None
    
    def __init__(self,name,abbrev,reference,endmember_names,\
                      formula_str,conversion_strs,\
                      **kwargs):
        # van Laar weights (nc)
        d     = kwargs.get('d', [])
        # interaction weights (nc x nc)
        W     = kwargs.get('W', [])
        sites = kwargs.get('sites', [])
        TC0   = kwargs.get('TC0', None)
        VD    = kwargs.get('VD', None)
        SD    = kwargs.get('SD', None)
        R     = kwargs.get('R', 8.31446261815324)
        T_r   = kwargs.get('T_r', 300.0)
        W_V   = kwargs.get('W_V', [])
        self.name   = name
        self.abbrev = abbrev
        self.reference = reference
        
        # endmembers (nc)
        self.endmember_names = endmember_names
        self.formula_str     = formula_str
        self.conversion_strs = conversion_strs

        # number of components
        self.nc = len(endmember_names)
        # numbers of sites
        self.ns = len(sites)
        # number of pairs
        self.np = self.nc*(self.nc-1)//2
        
        self.test_strs = ['['+repr(i)+'] > 0.0' for i in  range(self.nc)]

        # site multiplicities (ns)
        self.m = [sites[i][0] for i in range(self.ns)]
        # species coefficients for sites (ns x nc)
        self.C = [sites[i][1] for i in range(self.ns)]
        
        self.Rsym  = sym.symbols('R', real=True)
        self.Trsym = sym.symbols('T_r', real=True)
        self.dsyms = [sym.symbols('d_{}'.format(i), real=True) for i in range(len(d))]
        Wname = []
        Wsyms = []
        W_Vname = []
        W_Vsyms = []
        self.Wsarr = []
        if len(W) > 0:
            Wname = ['W_{}{}'.format(i,j) for j in range(2,self.nc+1) for i in range(1,j)]
            Wsyms = [sym.symbols('W_{}{}'.format(i,j), real=True) for j in range(2,self.nc+1) for i in range(1,j)]
            self.Wsarr = [[0 for j in range(self.nc)] for i in range(self.nc)]
            k = 0
            for i in range(1,self.nc):
                for j in range(i):
                    self.Wsarr[i][j] = Wsyms[k]
                    self.Wsarr[j][i] = Wsyms[k]
                    k = k + 1
            W_Vname = ['W_V_{}{}'.format(i,j) for j in range(2,self.nc+1) for i in range(1,j)]
            W_Vsyms = [sym.symbols('W_V_{}{}'.format(i,j), real=True) for j in range(2,self.nc+1) for i in range(1,j)]
            self.W_Vsarr = [[0 for j in range(self.nc)] for i in range(self.nc)]
            k = 0
            for i in range(1,self.nc):
                for j in range(i):
                    self.W_Vsarr[i][j] = W_Vsyms[k]
                    self.W_Vsarr[j][i] = W_Vsyms[k]
                    k = k + 1
        
        if TC0 is not None and not numpy.isnan(numpy.sum(TC0)):
            assert(VD is not None and not numpy.isnan(numpy.sum(VD)))
            assert(SD is not None and not numpy.isnan(numpy.sum(SD)))
            if not pd.api.types.is_list_like(TC0): TC0 = [TC0]
            if not pd.api.types.is_list_like(VD):  VD = [VD]
            if not pd.api.types.is_list_like(SD):  SD = [SD]
            assert(len(TC0)==self.nc)
            assert(len(VD)==self.nc)
            assert(len(SD)==self.nc)
            self.TC0 = TC0
            self.VD  = VD
            self.SD  = SD
            self.TC0syms = [sym.symbols('T_C0_{}'.format(i), real=True) for i in range(len(TC0))]
            self.VDsyms  = [sym.symbols('V_D_{}'.format(i), real=True)  for i in range(len(VD))]
            self.SDsyms  = [sym.symbols('S_D_{}'.format(i), real=True)  for i in range(len(SD))]
        
        self.param_strs = ['R', 'T_r']+\
                          ['d_{}'.format(i) for i in range(len(d))]+\
                          Wname
        self.param_unit = ['J/K/mol', 'K']+\
                          ['']*len(d)+\
                          ['J/mol']*len(Wsyms)
        self.param_syms = [self.Rsym, self.Trsym]+\
                          self.dsyms+\
                          W_Vsyms

        self.param_strs += W_Vname
        self.param_unit += ['J/mol/bar']*len(W_Vsyms)
        self.param_syms += W_Vsyms

        if self.TC0 is not None:
            self.param_strs += ['T_C0_{}'.format(i) for i in range(len(self.TC0))]
            self.param_strs += ['V_D_{}'.format(i) for i in range(len(self.VD))]
            self.param_strs += ['S_D_{}'.format(i) for i in range(len(self.SD))]
            self.param_unit += ['K']*len(self.TC0)
            self.param_unit += ['cm^3/mol/10']*len(self.VD)
            self.param_unit += ['J/mol/K']*len(self.SD)
            self.param_syms += self.TC0syms
            self.param_syms += self.VDsyms
            self.param_syms += self.SDsyms

        self.param_vals = {'R' : R, \
                           'T_r' : T_r}
        self.param_vals.update(dict([('d_{}'.format(i),d[i]) for i in range(len(d))]))
        self.param_vals.update(dict([(Wname[i],W[i]) for i in range(len(W))]))
        self.param_vals.update(dict([(W_Vname[i],W_V[i]) for i in range(len(W_V))]))
        if self.TC0 is not None:
            self.param_vals.update(dict([('T_C0_{}'.format(i),self.TC0[i]) for i in range(len(self.TC0))]))
            self.param_vals.update(dict([('V_D_{}'.format(i),self.VD[i]) for i in range(len(self.VD))]))
            self.param_vals.update(dict([('S_D_{}'.format(i),self.SD[i]) for i in range(len(self.SD))]))

        self.model = coder.SimpleSolnModel.from_type(nc=self.nc)
        
        self.n  = self.model.n
        self.nT = self.model.nT
        self.X  = self.n/self.nT
        
        self.T  = self.model.get_symbol_for_t()
        self.P  = self.model.get_symbol_for_p()
        self.mu = self.model.mu
    
    def G_ss_default(self):
        G_ss = (self.n.transpose()*self.mu)[0]
        return G_ss

    def G_config_default(self):
        if self.ns==0: return 0
        C = sym.Matrix(self.C)
        X_o = C*self.X
        M = sym.Matrix(self.m)
        G = self.nT*self.Rsym*self.T*M.T*sym.diag(*X_o)*X_o.applyfunc(sym.log)
        G_config = G[0].simplify()
        return G_config
    
    def G_excess_default(self):
        if len(self.Wsarr)==0: return 0
        W = sym.Matrix(self.Wsarr)
        for i in range(self.nc):
            for j in range(self.nc):
                W[i,j] = W[i,j]/(self.dsyms[i] + self.dsyms[j])
        Xdsum = sum([self.X[i]*self.dsyms[i] for i in range(self.nc)])
        Phi = sym.diag(*self.dsyms)*self.X/Xdsum
        G_excess = self.nT*Xdsum*Phi.T*W*Phi
        G_excess = G_excess[0].simplify()
        return G_excess
    
    def G_landau_default(self):
        if self.TC0 is None: return 0
        TC = [self.TC0syms[i] + self.VDsyms[i]*self.P/self.SDsyms[i] for i in range(len(self.TC0))]
        Q2 = [sym.sqrt(1-self.T/TC[i]) for i in range(len(TC))]
        G_landau_hip = [self.n[i]*self.SDsyms[i]*((self.T-TC[i])*Q2[i] + self.TC0syms[i]*Q2[i]*Q2[i]*Q2[i]/3) for i in range(len(TC))]
        G_landau_hip = [G_landau_hip[i].simplify() for i in range(len(TC))]
        G_landau = sum([sym.Piecewise((0, self.T >= TC[i]), (G_landau_hip[i], True)) for i in range(len(TC))])
        return G_landau
    
    def G_default(self):
        return self.G_ss_default() + self.G_config_default() + self.G_excess_default() + self.G_landau_default()
    def G_landau_2021(self):
        if self.TC0 is None: return 0

        """See SLB 2021 GJI.
        Three critical terms: TC0, VD, SD.
        These are zero for most endmembers, which can lead to nans and infinities if not careful.
        For species with Landau transitions:
            TC0 > 0 (commonly 5 K).
            SD > 0 (though vanishingly small for staurolite).
            VD >= 0 (but most commonly zero).
        """
        TC = [sym.S.Zero if self.TC0[i]==0 else self.TC0syms[i] + self.VDsyms[i]*self.P/self.SDsyms[i] for i in range(len(self.TC0))]


        # Following SLB 2021 appendix, Q has a maximum of 2 (Q**2  <= 4).
        Q2 = [sym.Piecewise((sym.S.Zero, self.T>=TC[i]), (sym.Min(sym.S(4), sym.sqrt((TC[i]-self.T)/self.TC0syms[i])), True)) for i in range(len(TC))]

        G_landau_hip = [sym.S.Zero if self.TC0[i]==0 else self.n[i]*self.SDsyms[i]*((self.T-TC[i])*(Q2[i]-sym.S(1)) + self.TC0syms[i]*(Q2[i]*Q2[i]*Q2[i]-sym.S(1))/sym.S(3)) for i in range(len(TC))]
        G_landau_hip = [G_landau_hip[i].simplify() for i in range(len(TC))]
        G_landau = sum(G_landau_hip)
        return G_landau
    
    def G_excess_2021(self):
        if len(self.Wsarr)==0: return 0
        W = sym.Matrix(self.Wsarr)
        W_V = sym.Matrix(self.W_Vsarr)
        for i in range(self.nc):
            for j in range(self.nc):
                # Following SLB2021, W_V is multiplied by pressure
                W[i,j] = (W[i,j] + W_V[i,j]*self.P)/(self.dsyms[i] + self.dsyms[j])
        Xdsum = sum([self.X[i]*self.dsyms[i] for i in range(self.nc)])
        Phi = sym.diag(*self.dsyms)*self.X/Xdsum
        G_excess = self.nT*Xdsum*Phi.T*W*Phi
        G_excess = G_excess[0].simplify()
        return G_excess
    
    def G_2021(self):
        return self.G_ss_default() + self.G_config_default() + self.G_excess_2021() + self.G_landau_2021()

    def add_parameter(self, name, unit, symbol, value):
        self.param_strs.append(name)
        self.param_unit.append(unit)
        self.param_syms.append(symbol)
        self.param_vals[name] = value
    
    def params(self):
        return list(zip(self.param_strs, self.param_unit, self.param_syms))

    def values_dict(self):
        values_dict = dict(
                           formula_string=self.formula_str,
                           conversion_string=self.conversion_strs,
                           test_string=self.test_strs,
                           name=self.name,
                           abbrev=self.abbrev,
                           reference=self.reference,
                           endmembers = self.endmember_names
                          )
        values_dict.update(self.param_vals)
        return values_dict
    
    def add_potential_to_model(self):
        self.model.add_potential_to_model('G', self.G, self.params())

    def set_model_values(self):
        values_dict = self.model.get_values()
        values_dict.update(self.values_dict())
        self.model.set_values(values_dict)
    
    def tofile(self,outdir):
        if self.G is None: self.G = self.G_default()
        self.add_potential_to_model()
        self.set_model_values()
        filename = self.model.to_xml(path=outdir)
        return filename

