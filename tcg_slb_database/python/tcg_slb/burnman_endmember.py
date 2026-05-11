import sympy as sym
import re
import types
from thermocodegen.coder import coder
import burnman

def get_names_for_database(database):
    db = getattr(burnman.minerals, database)
    potential_minerals = dir(db)
    name_dict = {}
    for m in potential_minerals:
        try:
            get_mineral = getattr(db, m)
            mineral = get_mineral()
            # check if it's already been set
            if mineral.params['name'] in name_dict:
                print(f"Warning: duplicate name {mineral.params['name']} for minerals {name_dict[mineral.params['name']]} and {m}. Skipping {m}.")
                # if this is the abbreviated name, skip it
                # but if this is the longer name, overwrite the dict
                if len(m) > len(name_dict[mineral.params['name']]):
                    name_dict[mineral.params['name']] = m
                else:
                    continue
            else:
                name_dict[mineral.params['name']] = m
        except Exception as e:
            print(f"Skipping {m} due to error: {e}")  
    names = list(name_dict.values())
    print(names)
    print(len(names))
    return list(name_dict.values())

class BurnmanEndmember:
    name    = None
    formula = None
    model   = None
    param_vals = None
    A       = None
    
    def __init__(self,database,name,reference,**kwargs):
        self.database = database
        self.name = name
        self.reference = reference
        self.model = coder.StdStateModel.from_type('TV')
        self.T = self.model.get_symbol_for_t()
        self.V = self.model.get_symbol_for_v()
        self.T_r = self.model.get_symbol_for_tr()
        self.V_r = self.model.get_symbol_for_vr()

        self.param_strs = [ ]
        self.param_unit = [ ]
        symdict = dict( ( p, sym.symbols(p,real=True) ) for p in self.param_strs )
        self.syms = types.SimpleNamespace(**symdict)
        self.param_syms = [ symdict[p] for p in self.param_strs ]
        required_params = self.param_strs

        missing_params = [ p for p in required_params if p not in kwargs ]
        if len(missing_params) > 0:
            raise Exception("Not all parameter values provided.  Missing: "+", ".join(missing_params))
        
        self.param_vals = dict( (p, kwargs[p]) for p in required_params )

        name = self.name.replace(f'_{self.database}_em','')
        db = getattr(burnman.minerals, self.database)
        get_mineral = getattr(db, name)
        self.mineral = get_mineral()


    def A_default(self):
        A = self.mineral.method.helmholtz_free_energy(1, self.T, self.V, self.mineral.params)
        return A

    def params(self):
        return list(zip(self.param_strs, self.param_unit, self.param_syms))

    def values_dict(self):
        print(self.mineral.params)
        values_dict = dict(
                           name=self.name,
                           formula=f"{''.join([f'{k}({v})' for k,v in self.mineral.formula.items()])}",
                           reference=self.reference
                          )
        values_dict.update(self.param_vals)
        return values_dict
    
    def set_model_values(self):
        values_dict = self.model.get_values()
        values_dict.update(self.values_dict())
        self.model.set_values(values_dict)
        print(values_dict)

    def add_potential_to_model(self):
        self.model.add_potential_to_model('A', self.A, self.params())
    
    def tofile(self,outdir):
        if self.A is None: self.A = self.A_default()
        self.add_potential_to_model()
        self.set_model_values()
        filename = self.model.to_xml(path=outdir,filename=self.name+".emml",verbose=True)
        return filename

