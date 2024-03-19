import sys, os

sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir, 'tcg_slb','python'))

from python.tcg import get_reaction,latex_reactions,get_names,x2c,phi2F,custom_solve
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.pyplot as plt
from pathlib import Path
from tcg_slb.base import GPa2Bar
from tcg_slb.phasediagram.scipy import ScipyPDReactiveODE
import multiprocessing as mp
from multiprocessing import Pool
from python.perplex import ppx_point_composition, ppx_profile_data

### ------------ INPUTS -------------------
reference= 'parallel_profile'
composition = 'hacker_2015_md_xenolith'
rxn_name = 'eclogitization_2024_stx21_rx'

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

# ============= Parse arguments for CLI =============

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--composition")
    parser.add_argument("-e", "--end_t")
    parser.add_argument("-r", "--rxn_name")
    parser.add_argument("-d", "--da")
    parser.add_argument("-n", "--num_procs")

    args = parser.parse_args()

    if args.composition is not None:
        composition = args.composition
    if args.end_t is not None:
        end_t = float(args.end_t)
    if args.rxn_name is not None:
        rxn_name = args.rxn_name
    if args.da is not None:
        Da = float(args.da)
    if args.num_procs is not None:
        processes = int(args.num_procs)

#====================================================

outputPath = Path("output",reference,composition,rxn_name)
outputPath.mkdir(parents=True, exist_ok=True)

df = ppx_profile_data(composition)
T_range = df['T(K)'].to_numpy() # K
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


rxn = get_reaction(rxn_name)
rxn.set_parameter("T0", Tr)
table = latex_reactions(rxn)

with open(Path(outputPath,"reaction_table_body.tex"), "w") as text_file:
    text_file.write(table)

Fi0, Xik0, phii0, cik0 = ppx_point_composition(rxn, composition)

phase_names, endmember_names = get_names(rxn)

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
    custom_solve(ode,T,GPa2Bar(P),Fi0,cik0,end_t,Da=Da,eps=eps,max_steps=max_steps)

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

fig = plt.figure(figsize=(10,10))

plt.plot(P_range,rho_final)
plt.suptitle(composition.replace("_", " ").capitalize())
plt.savefig(Path(outputPath,'density.png'))

fig = plt.figure(figsize=(10,5))
ax = plt.gca()
ax3 = ax.twiny()
df = ppx_profile_data(composition)

# T(K), P(bar), Pl, Pl, Cpx, Opx, qtz, Gt, ky  
phase_name_to_col_name = {
    "Clinopyroxene_stx21_ph":"Cpx3",
    "Orthopyroxene_stx21_ph":"Opx",
    "Quartz_stx21_ph":"qtz",
    "Feldspar_stx21_ph":"Pl2",
    "Garnet_stx21_ph":"Gt2",
    "Kyanite_stx21_ph":"ky",
    "Spinel_stx21_ph":"Sp",
    "Olivine_stx21_ph": "O"
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
if(df is not None):
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
            ax.legend(phase_names + ["(EQ)"])


if(df is not None):
    if("O" in df.columns and "Olivine_stx21_ph" not in [p.name() for p in rxn.phases()]):
        y = df["O"]/100
        ax.plot(xvar,y,"-",linewidth=1,alpha=0.5,color="black")
    if("Aki" in df.columns):
        y = df["Aki"]/100
        ax.plot(xvar,y,"-",linewidth=1,alpha=0.5,color="black")

#plt.suptitle(composition.replace("_", " ").capitalize())
plt.savefig(Path(outputPath,"phases.png"))

ax2 = ax.twinx()
ax2.set_ylabel('abs. diff.')  # we already handled the x-label with ax1
ax2.set_ylim([0.02,0])
if(df is not None):
    for i, phase in enumerate(rxn.phases()):
        pname = phase.name()
        col = phase_name_to_col_name[pname]
        if col not in df.columns:
            continue
        h = hs[i]
        y = df[col]/100
        y1 = phii_final[:,i]
        diff = y-y1
        ax2.plot(xvar, np.abs(diff),"-",color=h[-1].get_color(),linewidth=0.5)

plt.suptitle(composition.replace("_", " ").capitalize())
plt.savefig(Path(outputPath,"phases_errors.png"))           

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

for k in range(K):
    name = endmember_names[k]
    line_style = line_style_by_endmember[name]
    plt.plot(xvar,cik_final[:,k], line_style)

plt.ylim([0,1])
plt.xlim(xlimits)
plt.legend(endmember_names)
plt.suptitle(composition.replace("_", " ").capitalize())
plt.savefig(Path(outputPath,'endmembers_cik.png'))

fig = plt.figure(figsize=(12,12))
axi = fig.add_subplot(1,1,1)

for k in range(K):
    name = endmember_names[k]
    line_style = line_style_by_endmember[name]
    plt.plot(xvar,Xik_final[:,k], line_style)

plt.ylim([0,1])
plt.xlim(xlimits)
plt.legend(endmember_names)
plt.suptitle(composition.replace("_", " ").capitalize())
plt.savefig(Path(outputPath,'endmembers_Xik.png'))

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

plt.ylim([0,1])
plt.xlim(xlimits)
plt.legend(endmember_names)
plt.suptitle(composition.replace("_", " ").capitalize())
plt.savefig(Path(outputPath,'endmembers_wtpc.png'))