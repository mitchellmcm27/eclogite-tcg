import sys, os, time
import numpy as np
import sympy as sym
from scipy.integrate import solve_ivp
from scipy.integrate import solve_ivp, BDF, OdeSolution
from scipy.optimize import OptimizeResult as OdeResult
from tcg_slb.phasediagram.scipy import ScipyPDReactiveODE
from typing import TypedDict,List,Tuple

MESSAGES = {
    0: "The solver successfully reached the end of the integration interval.",
    1: "A termination event occurred.",
    2: "The solver failed to converge."
}

sym.init_printing()

def x2c(rxn, Xik0):
    return np.asarray([c for (i, ph) in enumerate(rxn.phases()) for c in ph.x_to_c(Xik0[i])])

def phi2F(rxn, phii, cik, T=900.,p=10000.,eps=1e-5):
    densities = []
    C = rxn.zero_C()
    Ki = 0
    for i,ph in enumerate(rxn.phases()):
        n = len(ph.endmembers())
        C[i] = cik[Ki:Ki+n]
        Ki = Ki+n
    # regularize C
    C = [np.maximum(np.asarray(C[i]), eps*np.ones(len(C[i]))) for i in range(len(C))]
    C = [np.asarray(C[i])/sum(C[i]) for i in range(len(C))]

    densities = [ph.rho(T, p, C[i]) for i,ph in enumerate(rxn.phases())]
    mass_tot = np.sum(np.asarray(densities) * np.asarray(phii))
    Fi = np.asarray([v*densities[i]/mass_tot for (i, v) in enumerate(phii)])
    return Fi

def get_reaction(rxnName):
    pv = repr(sys.version_info.major)+'.'+repr(sys.version_info.minor)
    path = os.path.join(os.path.pardir, 'tcg_slb_database','database', 'install', rxnName,
                        'lib', 'python'+pv, 'site-packages/')  # the final slash is necessary
    sys.path.append(path)
    tcgdb = __import__('py_'+rxnName)
    importer = getattr(tcgdb, rxnName)
    rxn = importer()
    return rxn

def get_names(rxn):
    phase_names = [p.name() for p in rxn.phases()]
    phase_names = [s.replace("_stx21_ph","") for s in phase_names]

    endmember_names = [em.name() for p in rxn.phases() for em in p.endmembers()]
    endmember_names = [s.replace("_stx21_em","") for s in endmember_names]
    return phase_names, endmember_names

def latex_reactions(rxn):
    '''
    rudnick_2014_lower_crust
    Reaction object: eclogitization_agu18_stx21_rx

    Phase 0 Clinopyroxene_stx21_ph (cpx)
         Endmember 0 Diopside_stx21_em : CaMgSi2O6_(cpx)
         Endmember 1 Hedenbergite_stx21_em : CaFeSi2O6_(cpx)
         Endmember 2 Clinoenstatite_stx21_em : Mg2Si2O6_(cpx)
         Endmember 3 CaTschermaks_stx21_em : CaAl2SiO6_(cpx)
         Endmember 4 Jadeite_stx21_em : NaAlSi2O6_(cpx)
    Phase 1 Orthopyroxene_stx21_ph (opx)
         Endmember 0 Enstatite_stx21_em : Mg2Si2O6_(opx)
         Endmember 1 Ferrosilite_stx21_em : Fe2Si2O6_(opx)
         Endmember 2 MgTschermaks_stx21_em : MgAl2SiO6_(opx)
         Endmember 3 OrthoDiopside_stx21_em : CaMgSi2O6_(opx)
    Phase 2 Quartz_stx21_ph (qtz)
         Endmember 0 Quartz_stx21_em : SiO2_(qtz)
    Phase 3 Feldspar_stx21_ph (plg)
         Endmember 0 Anorthite_stx21_em : CaAl2Si2O8_(plg)
         Endmember 1 Albite_stx21_em : NaAlSi3O8_(plg)
    Phase 4 Garnet_stx21_ph (gt)
         Endmember 0 Pyrope_stx21_em : Mg3Al2Si3O12_(gt)
         Endmember 1 Almandine_stx21_em : Fe3Al2Si3O12_(gt)
         Endmember 2 Grossular_stx21_em : Ca3Al2Si3O12_(gt)
         Endmember 3 MgMajorite_stx21_em : Mg4Si4O12_(gt)
         Endmember 4 NaMajorite_stx21_em : Na2Al2Si4O12_(gt)
    Phase 5 Kyanite_stx21_ph (ky)
         Endmember 0 Kyanite_stx21_em : Al2SiO5_(ky)

    Reaction 0
         0.666667 CaFeSi2O6_(cpx) + 0.333333 Mg2Si2O6_(opx) -> 0.666667 CaMgSi2O6_(cpx) + 0.333333 Fe2Si2O6_(opx)
    Reaction 1
         0.6 Fe2Si2O6_(opx) + 0.4 Mg3Al2Si3O12_(gt) -> 0.6 Mg2Si2O6_(opx) + 0.4 Fe3Al2Si3O12_(gt)
    Reaction 2
         0.75 CaFeSi2O6_(cpx) + 0.25 Mg3Al2Si3O12_(gt) -> 0.75 CaMgSi2O6_(cpx) + 0.25 Fe3Al2Si3O12_(gt)
    Reaction 3
         Mg2Si2O6_(opx) -> Mg2Si2O6_(cpx)
    '''
    
    import io
    from contextlib import redirect_stdout

    with io.StringIO() as buf, redirect_stdout(buf):
        rxn.report()
        report = buf.getvalue()

    lines = report.splitlines()
    ref = lines[0]
    name = lines[1]
    phases=[]
    reactions=[]
    names = []
    for i in range(2,len(lines)):
        line = lines[i]
        if(line[0:5]=="Phase"):
            phases.append(line[8:].strip())
        if(line[0:8]=="Reaction"):
            n = int(line[8:].strip())+1
            names.append("\\textbf{"+str(n)+".}")
            rxn = lines[i+1].strip()
            rxn = rxn.replace("0.333333","1/3")
            rxn = rxn.replace("0.666667","2/3")
            rxn = rxn.replace("0.75","3/4")
            rxn = rxn.replace("0.5","1/2")
            rxn = rxn.replace("0.25","1/4")
            rxn = rxn.replace("0.2","1/5")
            rxn = rxn.replace("0.4","2/5")
            rxn = rxn.replace("0.6","3/5")
            rxn = rxn.replace("0.8","4/5")
            rxn = rxn.replace("0.166667","1/6")
            rxn = rxn.replace("_","")
            rxn = rxn.replace("->", "=")
            rxn = "\ce{" + rxn + "}"
            reactions.append(rxn)

    # 1 & $\ce{2/3 CaFeSi2O6}^\text{cpx} + \ce{1/3 Mg2Si2O6}^\text{opx} = \ce{2/3 CaMgSi2O6}^\text{cpx} + \ce{1/3 Fe2Si2O6}^\text{opx}$ \\
    table = "\n".join(["{} & {} \\\\".format(names[i], reactions[i]) for i in range(len(names))])
    return table

def composition_to_label(c):
    default = c.replace("_", " ").capitalize()
    dic = {
        "hacker_2015_md_xenolith": "Hacker et al. (2015) median xenolith",
        "hacker_2015_bin_1": "Hacker et al. (2015) bin-1",
        "hacker_2015_bin_2": "Hacker et al. (2015) bin-2",
        "hacker_2015_bin_3": "Hacker et al. (2015) bin-3",
        "hacker_2015_bin_4": "Hacker et al. (2015) bin-4",
        "sammon_2021_lower_crust": "Sammon & McDonough (2021) lower crust",
        "sammon_2021_deep_crust": "Sammon & McDonough (2021) deep crust",
        "xu_2008_basalt": "Xu et al. (2008) basalt",
        "xu_2008_pyrolite": "Xu et al. (2008) pyrolite",
        "xu_2008_harzburgite": "Xu et al. (2008) harzburgite",
        "zhang_2022_cd07-2": "Zhang et al. (2022) granulite",
        "zhang_2006_mafic_granulite": "Zhang et al. (2006) granulite",
        "bhowany_2018_hol2a": "Bhowany et al. (2018) eclogite (hol2a)"
    }
    return dic.get(c, default)

def custom_solve(ode:ScipyPDReactiveODE, T:float, p:float, mi0:List[float], Cik0:List[float], end:float, **kwargs):

    max_steps = kwargs.pop("max_steps", np.inf)

    ode.set_initial_params(T,p,mi0,Cik0,**kwargs)

    rtol   = kwargs.get('rtol', 1.e-5)
    atol   = kwargs.get('atol', 1.e-9)

    u0 = np.concatenate((mi0,Cik0))

    tic = time.perf_counter()

    t0 = 0.
    solver = BDF(ode.rhs, t0, u0, float(end), rtol=rtol, atol=atol, jac=ode.jac)
    ts = [t0]
    ys = [u0]
    interpolants = []
    status = None
    num_steps = 0
    while status is None:
        message = solver.step()
    
        if solver.status == 'finished':
            status = 0
        elif solver.status == 'failed':
            status = -1
            break

        if num_steps == max_steps:
            status = 2
        
        t_old = solver.t_old
        t = solver.t
        y = solver.y
        output = solver.dense_output()
        interpolants.append(output)
        ts.append(t)
        ys.append(y)
        num_steps = num_steps + 1
        
    message = MESSAGES.get(status, message)
    ts = np.array(ts)
    ys = np.vstack(ys).T
    ode_solution = OdeSolution(ts, interpolants)
    sol = OdeResult(t=ts, y=ys, sol=ode_solution,
            nfev=solver.nfev, njev=solver.njev, nlu=solver.nlu,
            status=status, message=message, success=status >= 0)
    ode.sol = sol
    toc = time.perf_counter()
    ode.stime = toc-tic