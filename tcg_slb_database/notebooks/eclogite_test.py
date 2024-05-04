import sys, os
import multiprocessing as mp
from multiprocessing import Pool
import pandas as pd
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir,'python'))

import numpy as np
from matplotlib import pyplot as plt
from pathlib import Path
from tcg_slb.base import GPa2Bar
from tcg_slb.phasediagram.scipy import ScipyPDReactiveODE
from scipy.integrate import solve_ivp
from io import StringIO

### ------------ INPUTS -------------------
reference= 'eclogite_test'
rxn_name = 'eclogitization_2024_slb21_rx'

# end time of reactions, change with -e argument
end_t = 1
# reaction's characteristic temperature (T_r)
Tr = 5500.+273.15 

# only phases greater than this fraction will be plotted
phasetol = 1.e-3 # 1.e-2

# Damkhoeler number, change with -d argument
Da = 1e6 # 1.0
# regularization parameter for compositions
eps = 1.e-5 # 1.e-2
# these numbers seem to work very well with eps = 1e-5??
rtol = 1.e-5 # relative tolerance, default 1e-5
atol = 1.e-9 # absolute tolerance, default 1e-9

# large number
max_steps = 4e3 # 4e3 is reasonable

# number of processes, edit with -n argument
processes = mp.cpu_count()

# ------------------------------------------

outputPath = Path("figures",reference)
outputPath.mkdir(parents=True, exist_ok=True)

pv = repr(sys.version_info.major)+'.'+repr(sys.version_info.minor)
path = os.path.join(os.path.pardir, 'database', 'install', rxn_name,
                    'lib', 'python'+pv, 'site-packages/')  # the final slash is necessary
print(path)
sys.path.append(path)
tcgdb = __import__('py_'+rxn_name)
importer = getattr(tcgdb, rxn_name)
rxn = importer()

def ppx_point_composition(rxn):
    filepath = "eclogite_eq_point.txt"
    print(filepath)
    try:
        x = "" # mol fractions
        F = "" # mass fractions
        with open(filepath) as file:
            copy = False
            for line in file:
                if line.strip().startswith("Phase Compositions"):
                    copy = True
                elif line.strip().startswith("Phase speciation"):
                    copy = False
                elif copy:
                    F += line.strip() + "\n"
        with open(filepath) as file:
            copy = False
            for line in file:
                if line.strip().startswith("Phase speciation"):
                    copy = True
                elif line.strip().startswith("Structural formulae for 688 format solution models"):
                    copy = False
                elif copy:
                    x += line.strip() + "\n"
  
        F = F.strip().split("\n")
        F[0] = F[0].replace(" %", "%")
        F = "\n".join(F)
        F_df = pd.read_csv(StringIO(F), delimiter="\s+", header=0)
        Fs = F_df['wt%']
        Fi0 = [0 for p in rxn.phases()]

        for n,p in enumerate(rxn.phases()):
            if p.abbrev()=="cpx":
                Fi0[n] = Fs.get("Cpx",0.)
            if p.abbrev()=="opx":
                Fi0[n] = Fs.get("Opx",0.)
            if p.abbrev()=="qtz":
                Fi0[n] = Fs.get("qtz",0.)
            if p.abbrev()=="plg":
                Fi0[n] = Fs.get("Pl",0.)
            if p.abbrev()=="gt":
                Fi0[n] = Fs.get("Gt",0.)
            if p.abbrev()=="ky":
                Fi0[n] = Fs.get("ky",0.)
            if p.abbrev()=="sp":
                Fi0[6] = Fs.get("Sp",0.)
            if p.abbrev()=="co":
                Fi0[6] = Fs.get("Aki",0.)

        Fi0 = [f/100 for f in Fi0]
        
        x = x.strip().split("\n")
        Xik0 = rxn.zero_C()
        for i,c in enumerate(Xik0):
            Xik0[i][0] = 1.

        import re
        for line in x:
            l = re.sub(r'\s+',' ', line.strip())
            [phase, rest] = l.split(" ", 1)
            ems = []
            for vals in rest.split(','):
                [k,v] = vals.strip().split(" ")
                ems.append(float(v))
            if(phase=="Pl"):
                Xik0[3] = [ems[1],ems[0]]
            elif(phase=="Sp"):
                if len(Xik0)>=7:
                    Xik0[6] = [ems[1],ems[0]]
            elif(phase=="Cpx"):
                Xik0[0] = [ems[1], ems[2], ems[3], ems[4], ems[0]]
            elif(phase=="Opx"):
                Xik0[1] = [ems[1], ems[2], ems[3], ems[0]]
            elif(phase=="Gt"):
                if(len(ems)==5):
                    # 5 endmember garnet
                    Xik0[4] = [ems[4], ems[3], ems[2], ems[1], ems[0]]
                    #print(Xik0)
                    # regularize 3-component garnet
                    g3 = (1-(Xik0[4][0]+Xik0[4][1]+Xik0[4][2]))/3.0
                    Xik0[4][0] += g3
                    Xik0[4][1] += g3
                    Xik0[4][2] += g3
                    Xik0[4][3] = 0.0
                    Xik0[4][4] = 0.0

                elif(len(ems)==3):
                    # 3 endmember garnet
                    Xik0[4] = [ems[2], ems[1], ems[0], 0., 0.,]

        phii0 = None
        Cik0 = None
        return Fi0, Xik0, phii0, Cik0
    except Exception as e:
        print("There was an error:")
        print(e)
        rase(e)

def ppx_profile_data():
    filepath = "eclogite_eq_profile.dat"
    try:
        df = pd.read_csv(filepath, delimiter='\s+', skiprows=8, header=0)
        print(df)
        df.fillna(0, inplace=True)
        pl0 = df["Pl"]
        pl1 = 0 if not "Pl.1" in df.columns else df["Pl.1"]
        cpx0 = df["Cpx"]
        cpx1 = 0 if not "Cpx.1" in df.columns else df["Cpx.1"]
        cpx2 = 0 if not "Cpx.2" in df.columns else df["Cpx.2"]
        gt0 = df["Gt"]
        gt1 = 0 if not "Gt.1" in df.columns else df["Gt.1"]
        df["Pl2"] = pl0 + pl1
        df["Cpx3"] = cpx0 + cpx1 + cpx2
        df["Gt2"] = gt0 + gt1
        return df
    except Exception as e:
        print("Error parsing perple_x profile:")
        print(e)
        raise(e)

df = ppx_profile_data()

T_range = df['T(K)'].to_numpy()
P_range = df['P(bar)'].to_numpy()/1e4 # GPa
Tmin = np.min(T_range)
Tmax = np.max(T_range)
Pmin = np.min(P_range)
Pmax = np.max(P_range)

tdiff = np.abs(Tmax-Tmin)/500.
pdiff = np.abs(Pmax-Pmin)/2.

if(pdiff > tdiff):
    xaxis = "pressure"
    xvar = P_range
    xlabel = "Pressure (GPa)"
    xlimits = [Pmin, Pmax]

    x2var = T_range-273.15
    x2label = "Temperature (°C)"
    x2limits = [Tmin-273.15, Tmax-273.15]
else:
    xvar = T_range-273.15
    xaxis = "temperature"
    xlabel = "Temperature (°C)"
    xlimits = [Tmin-273.15, Tmax-273.15]

    x2var = P_range
    x2label = "Pressure (GPa)"
    x2limits = [Pmin, Pmax]

rxn.set_parameter("T0", Tr)

Fi0, Xik0, phii0, cik0 = ppx_point_composition(rxn)
phase_names = [p.name() for p in rxn.phases()]
phase_names = [s.replace("_slb21_ph","") for s in phase_names]
endmember_names = [em.name() for p in rxn.phases() for em in p.endmembers()]
endmember_names = [s.replace("_slb21_em","") for s in endmember_names]

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

cik0 = x2c(rxn, Xik0) if cik0 is None else cik0
Fi0 = phi2F(rxn, phii0, cik0) if Fi0 is None else Fi0

print([p.abbrev() for p in rxn.phases()])
print(cik0)
print(phase_names)

I = len(rxn.phases())
Kis = [len(rxn.phases()[i].endmembers()) for i in range(I)] # list, num EMs in each phase
K = sum(Kis)

rho_final = np.zeros(T_range.shape) # overall density
phases_final = ['' for i in T_range] # names of phases present
phii_final = np.empty(T_range.shape+(I,)) # phase vol. fractions
cik_final = np.empty(T_range.shape+(K,)) # EM mass fractions
Xik_final = np.empty(cik_final.shape) # EM mol fractions

# function to run in parallel
def task(i):
    P = P_range[i]
    T = T_range[i]

    ode = ScipyPDReactiveODE(rxn)
    #custom_solve(ode,T,GPa2Bar(P),Fi0,cik0,end_t,Da=Da,eps=eps,max_steps=max_steps)
    ode.solve(T,GPa2Bar(P),Fi0,cik0,end_t,Da=Da,eps=eps)
    rho = ode.final_rho() + 0.3

    cik = ode.sol.y[ode.I:ode.I+ode.K,-1]
    Fi = ode.sol.y[:ode.I,-1] # -1 = final time step
    C = ode.reshapeC(cik)
    Cs = ode.regularizeC(C)
    cik_reg = [c for arr in Cs for c in arr]
    rhoi = rxn.rho(T, GPa2Bar(P), Cs)
    v = Fi/rhoi
    vi  = 1./v.sum()
    phii = vi*Fi/rhoi
    Xik = [x for xarr in ode.rxn.C_to_X(C) for x in xarr]
    odephasenames, phaseabbrev = ode.final_phases(phasetol)
    phases = '+'.join(phaseabbrev)
    return rho, phases, phii, Fi, cik_reg, Xik, i

# Map blocks to block-level calculations
with Pool(processes,maxtasksperchild=12) as pool:
    # blocks until all finished
    sols = pool.map(task, range(len(T_range)))

for rho, phases, phii, Fi, cik, Xik, i in sols:
    rho_final[i] = rho
    phases_final[i] = phases
    phii_final[i,:] = phii
    cik_final[i,:] = cik
    Xik_final[i,:] = Xik

fig = plt.figure(figsize=(10,5))
ax = plt.gca()
ax3 = ax.twiny()

phase_name_to_col_name = {
    "Clinopyroxene_slb21_ph":"Cpx3",
    "Orthopyroxene_slb21_ph":"Opx",
    "Quartz_slb21_ph":"qtz",
    "Feldspar_slb21_ph":"Pl2",
    "Garnet_slb21_ph":"Gt2",
    "Kyanite_slb21_ph":"ky",
    "Spinel_slb21_ph":"Sp",
    "Olivine_slb21_ph": "O"
}

hs = []

for i, phase in enumerate(rxn.phases()):
    h = ax.plot(xvar,phii_final[:,i],'-', linewidth=3,alpha=0.5)
    hs.append(h)

ax.set_xlabel(xlabel)
ax.set_xlim(xlimits)
ax3.set_xlabel(x2label)
ax3.set_xlim(x2limits)
ax.set_ylabel("Phase vol. mode")
#ax.set_ylim([0,1])

for i, phase in enumerate(rxn.phases()):
    pname = phase.name()
    col = phase_name_to_col_name[pname]
    if col not in df.columns:
        print(col)
        continue
    h = hs[i]
    y = df[col]/100
    y1 = phii_final[:,i]
    diff = y-y1
    ax.plot(xvar,y,"--",linewidth=1, color=h[-1].get_color())
    if(i==0):
        ax.legend(phase_names + ["(Perple_X)"])

if("O" in df.columns and "Olivine_slb21_ph" not in [p.name() for p in rxn.phases()]):
    y = df["O"]/100
    ax.plot(xvar,y,"-",linewidth=1,alpha=0.5,color="black")
if("Aki" in df.columns):
    y = df["Aki"]/100
    ax.plot(xvar,y,"-",linewidth=1,alpha=0.5,color="black")

plt.savefig(Path(outputPath,"phases.png"))

fig = plt.figure(figsize=(12,12))
axi = fig.add_subplot(1,1,1)

line_style_by_endmember = {
    "Quartz": "-",
    "Kyanite": "-",
    "Diopside": ":",
    "Hedenbergite": ":",
    "Jadeite": ":",
    "CaTschermaks": ":",
    "Clinoenstatite": ":",
    'Enstatite':"-.", 
    'Ferrosilite':"-.", 
    'MgTschermaks':"-.", 
    'OrthoDiopside':"-.",
    "Anorthite": "-",
    "Albite": "-",
    'Pyrope': "--", 
    'Almandine':"--", 
    'Grossular':"--", 
    'MgMajorite':"--", 
    'NaMajorite':"--",
    'MgSpinel':"-",
    'Hercynite':"-",
    "Forsterite":":",
    "Fayalite":":"
    }

fig = plt.figure(figsize=(12,12))
axi = fig.add_subplot(1,1,1)

N = 0
for i, phase in enumerate(rxn.phases()):
    phii = phii_final[:,i]
    for k, em in enumerate(phase.endmembers()):
        em_c = cik_final[:,N+k]
        name = endmember_names[N+k]
        line_style = line_style_by_endmember[name]
        plt.plot(xvar,em_c*phii, line_style)
    N = N + len(phase.endmembers())

#plt.ylim([0,1])
plt.xlim(xlimits)
plt.legend(endmember_names)
axi.set_ylabel("Endmember compositions (wt%)")
plt.savefig(Path(outputPath,'endmembers_wtpc.png'))
