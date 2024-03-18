import sympy as sym
import re
import types
from thermocodegen.coder import coder

class SLBEndmember:
    name    = None
    formula = None
    model   = None
    param_vals = None
    A       = None
    
    def __init__(self,name,formula,reference,**kwargs):
        self.name = name
        self.formula = formula
        self.reference = reference
        self.model = coder.StdStateModel.from_type('TV')
        self.T = self.model.get_symbol_for_t()
        self.V = self.model.get_symbol_for_v()
        self.T_r = self.model.get_symbol_for_tr()
        self.V_r = self.model.get_symbol_for_vr()
        
        self.param_strs = [ 'a0', 'n', 'v0', \
                            'k00', 'k0p', 'theta0', \
                            'gamma0', 'q', 'r_Fe', \
                            'R' ]
        self.param_unit = [ 'J/m', '', 'J/bar-m',\
                            'bar', '', 'K', \
                            '', '', '', \
                            'J/K-m' ]
        symdict = dict( ( p, sym.symbols(p,real=True) ) for p in self.param_strs )
        self.syms = types.SimpleNamespace(**symdict)

        self.param_syms = [ symdict[p] for p in self.param_strs ]
        
        required_params = self.param_strs+['T_r','V_r']
        required_params.remove('R')    # R isn't actually required
        required_params.remove('r_Fe') # neither is r_Fe
        
        missing_params = [ p for p in required_params if p not in kwargs ]
        if len(missing_params) > 0:
            raise Exception("Not all parameter values provided.  Missing: "+", ".join(missing_params))
        
        self.param_vals = dict( (p, kwargs[p]) for p in required_params )
        self.param_vals['r_Fe'] = self.r_Fe()
        self.param_vals['R'] = kwargs.get('R', 8.31446261815324)
    
    def r_Fe(self):
        """Return number of Fe atoms in a formula."""
        formula_split = re.findall('[A-Z][^A-Z]*', self.formula)
        elem_count = {}
        for elem in formula_split:
            elem_split = re.split('(\d+)', elem)
            elem = elem_split[0].rstrip('(')
            elem_num = 1 if len(elem_split) == 1 else int(elem_split[1])
            elem_count[elem] = elem_count.get(elem,0) + elem_num
        return elem_count.get('Fe',0)

    def A_iso_default(self):
        """Isotropic strain contributions."""
        f = ((self.syms.v0/self.V)**(sym.S(2)/sym.S(3)) - 1)/2
        c1 = sym.S(9)/sym.S(2)*self.syms.k00*self.syms.v0
        c2 = self.syms.k0p-4
        A_iso = self.syms.a0 + c1 * f**2 * (1 + c2*f)
        return A_iso
    
    def A_quasi_default(self):
        db = coder.get_external_functions()['Debye']
        f = ((self.syms.v0/self.V)**(sym.S(2)/sym.S(3)) - 1)/2
        c3 = 6 * self.syms.gamma0
        c4 = c3 * (-2 + c3 - 3*self.syms.q)/2
        d0 = self.syms.theta0 * sym.sqrt( 1 + c3*f + c4 * f**2) 
        x = d0/self.T
        A_db = self.syms.n*self.syms.R*self.T*(3*sym.log(1-sym.exp(-x)) - db(x))
        A_quasi = A_db - A_db.subs(self.T,self.T_r)
        return A_quasi
    
    def A_mag_default(self):
        A_mag = -self.syms.R*self.T*self.syms.r_Fe*sym.log(5)
        return A_mag
    
    def A_default(self):
        return self.A_iso_default() + self.A_quasi_default() + self.A_mag_default()

    def add_parameter(self, name, unit, symbol, value):
        self.param_strs.append(name)
        self.param_unit.append(unit)
        self.param_syms.append(symbol)
        self.param_vals['name'] = value
    
    def params(self):
        return list(zip(self.param_strs, self.param_unit, self.param_syms))

    def values_dict(self):
        values_dict = dict(
                           name=self.name,
                           formula=self.formula,
                           reference=self.reference
                          )
        values_dict.update(self.param_vals)
        return values_dict
    
    def add_potential_to_model(self):
        self.model.add_potential_to_model('A', self.A, self.params())

    def set_model_values(self):
        values_dict = self.model.get_values()
        values_dict.update(self.values_dict())
        self.model.set_values(values_dict)
    
    def tofile(self,outdir):
        if self.A is None: self.A = self.A_default()
        self.add_potential_to_model()
        self.set_model_values()
        filename = self.model.to_xml(path=outdir)
        return filename

