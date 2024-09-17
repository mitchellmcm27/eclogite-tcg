import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir, 'tcg_slb_database','python'))

import pickle
from python.tcg import get_reaction,get_names,x2c,phi2F
import numpy as np
import numpy.ma as ma
from matplotlib import pyplot as plt
import matplotlib.pyplot as plt
from pathlib import Path
from multiprocessing import Pool
import multiprocessing as mp
from python.perplex import model_pyrolite_rho_gcc, ppx_point_composition
from scipy.integrate import solve_ivp
from geotherm_steady import geotherm_steady
import csv
from typing import TypedDict,List,Tuple

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

class TectonicSetting(TypedDict):
    setting:str
    L0:float
    As:float
    hr0:float
    k:float
    Ts:float
    Tlab:float

class InputScenario(TectonicSetting):
    Da:float
    composition:str
    T0:float
    T1:float
    qs0:float
    qs1:float
    P0:float
    cik0:List[float]
    Fi0:List[float]
    rho0:float

class OutputScenario(InputScenario):
    T: List[float] # K
    P: List[float] # bar
    rho: List[float]
    Fi: List[List[float]]
    cik: List[List[float]]
    Xik: List[List[float]]
    z: List[float]
    time: List[float]

####################
# Unit conversions #
####################

yr = 3.154e7
kyr = 1e3*yr
Myr = 1e6*yr
s = 1
mm = 1e-3
km = 1e3
g = 1e-3
cm = 1e-2

##########
# Inputs #
##########

save_output = False
load_output = True

reference = "parallel_experiment2"
rxn_name = "eclogitization_2024_slb21_rx"

# only phases greater than this fraction will be plotted
phasetol = 1.e-5 # default 1.e-2

# regularization parameter for compositions
eps = 1.e-5 # default 1.e-2
# these seem to work well with eps = 1e-5
rtol = 1.e-5 # relative tolerance
atol = 1.e-9 # absolute tolerance

z0  = 30. * km # initial Moho depth
v0 = 1.0 * mm/yr # Moho descent rate, m/s
h0 = 50. * km # total crustal thickening
z1 = z0 + h0 # final Moho depth
t0 = h0 / v0 # seconds
Tr = 5500.+273.15 # reaction's characteristic temperature (T_r)
crustal_rho = 2780.
gravity = 9.81
Da_eq = 1e6

# multiprocessing
num_processes =  mp.cpu_count()

# allows deterministic PDFs
pdf_metadata = {'CreationDate': None}

# Damkoehler numbers
Das = [1e-2, 3e-2, 1e-1, 3e-1, 1e0, 3e0, 1e1, 3e1, 1e2, 3e2, 1e3, 3e3, 1e4, 3e4, 1e5, 3e5, 1e6]
print(Das)

# default end time (scaled) is 1
end_t = 1.

# prefix file path for saving plots
prefix = None

# Compositions
compositions = [
    "sammon_2021_lower_crust",
    "sammon_2021_deep_crust",
    "hacker_2015_md_xenolith",
    "mackwell_1998_maryland_diabase"
]

color_by_composition = {
    'mackwell_1998_maryland_diabase': '#00adee', #blue
    'hacker_2015_md_xenolith': '#3cb371', # green
    'sammon_2021_lower_crust': '#be1e2d', # red
    'sammon_2021_deep_crust': '#f6921e' # yellow
}

tectonic_settings: List[TectonicSetting] = [
    {
        "setting": "a",
        "L0": 55.e3,
        "As": 2.0e-6,
        "hr0": 13.e3,
        "k": 3.0,
        "Ts": 10. + 273.15,
        "Tlab": 1330. + 273.15
    },    
    {
        "setting": "b",
        "L0": 60.e3,
        "As": 1.95e-6,
        "hr0": 12.5e3,
        "k": 3.0,
        "Ts": 10. + 273.15,
        "Tlab": 1330. + 273.15
    },
    {
        "setting": "c",
        "L0": 65.5e3,
        "As": 1.9e-6,
        "hr0": 12.e3,
        "k": 3.0,
        "Ts": 10. + 273.15,
        "Tlab": 1330. + 273.15
    },
    {
        "setting": "d",
        "L0": 71.5e3,
        "As": 1.85e-6,
        "hr0": 11.5e3,
        "k": 3.0,
        "Ts": 10. + 273.15,
        "Tlab": 1330. + 273.15
    },
    {
        "setting": "e",
        "L0": 78.e3,
        "As": 1.8e-6,
        "hr0": 11.0e3,
        "k": 3.0,
        "Ts": 10. + 273.15,
        "Tlab": 1330. + 273.15
    },
    {
        "setting": "f",
        "L0": 85.e3,
        "As": 1.75e-6,
        "hr0": 10.5e3,
        "k": 3.0,
        "Ts": 10. + 273.15,
        "Tlab": 1330. + 273.15
    },
    {
        "setting": "g",
        "L0": 92.5e3,
        "As": 1.7e-6,
        "hr0": 10.e3,
        "k": 3.0,
        "Ts": 10. + 273.15,
        "Tlab": 1330. + 273.15
    },
    {
        "setting":"h",
        "L0": 102.5e3,
        "As": 1.65e-6,
        "hr0": 9.5e3,
        "k": 3.0,
        "Ts": 10. + 273.15,
        "Tlab": 1330. + 273.15
    },
    {
        "setting":"i",
        "L0": 111.e3,
        "As": 1.6e-6,
        "hr0": 9.e3,
        "k": 3.0,
        "Ts": 10. + 273.15,
        "Tlab": 1330. + 273.15
    }
]


#######################
# Parse CLI arguments #
#######################

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--composition")
    parser.add_argument("-r", "--rxn_name")
    parser.add_argument("-q", "--quick", default=False, action="store_true")
    parser.add_argument("-n", "--num_processes")
    parser.add_argument("-f", "--force", default=False, action="store_true")
    parser.add_argument("-e", "--end_time")
    parser.add_argument("-p", "--prefix")

    args = parser.parse_args()

    if args.end_time is not None:
        print("Using custom end time of {}".format(args.end_time))
        end_t = float(args.end_time)
    if args.composition is not None:
        print("Using composition {}".format(args.composition))
        compositions = [args.composition]
    if args.rxn_name is not None:
        print("Using reaction {}".format(args.rxn_name))
        rxn_name = args.rxn_name
    if args.quick:
        print("Running in quick mode")
        Das = [d for d in Das if d <= 1e4]
        print("Only running Da = {}".format(Das))
    if args.num_processes is not None:
        num_processes = int(args.num_processes)
    if args.force:
        save_output = True
        load_output = False
    if args.prefix is not None:
        prefix = args.prefix

##########################
# Ready output directory #
##########################

output_path = Path("output",reference,rxn_name,prefix) if prefix is not None else Path("output",reference,rxn_name)
output_path.mkdir(parents=True, exist_ok=True)
pickle_path = Path(output_path,"_outs.pickle")

######################################
# Plot pressure-temperature profiles #
######################################

fig = plt.figure(figsize=(5,7))
cmap1 = plt.cm.get_cmap("coolwarm_r")
ax1 = plt.gca()
ax1.set_prop_cycle(plt.cycler("color", cmap1(np.linspace(0., 1., len(tectonic_settings)))))

geotherm_latex_table = {
    "heading": ["" for s in tectonic_settings],
    "As": np.zeros(len(tectonic_settings)),
    "L0": np.zeros(len(tectonic_settings)),
    "z0": np.zeros(len(tectonic_settings)),
    "hr0": np.zeros(len(tectonic_settings)),
    "qs0": np.zeros(len(tectonic_settings)),
    "qs1": np.zeros(len(tectonic_settings)),
    "T0": np.zeros(len(tectonic_settings)),
    "T1": np.zeros(len(tectonic_settings)),
}

for num_setting, setting in enumerate(tectonic_settings):
    L0 = setting["L0"]
    shortening = 1
    hr0 = setting["hr0"]
    conductivity = setting["k"]
    Ts = setting["Ts"]
    Tlab = setting["Tlab"]
    As = setting["As"]

    depths = np.linspace(0,1)
    depths_sc = L0*depths
    P =  depths_sc * crustal_rho * gravity / 1e5
    T, qs0 = geotherm_steady(depths,
                        L0*shortening,
                        shortening,
                        Ts=Ts,
                        Tlab=Tlab,
                        k=conductivity,
                        A=As,
                        hr0=hr0)


    p = plt.plot(T-273.15, depths_sc/1e3,linewidth=1, alpha=0.5)
    color = plt.gca().lines[-1].get_color()
    
    T0, _qs = geotherm_steady(z0/L0,
                    L0*shortening,
                    shortening,
                    Ts=Ts,
                    Tlab=Tlab,
                    k=conductivity,
                    A=As,
                    hr0=hr0)
    plt.plot(T0-273.15, z0/1e3,'.',color=color,alpha=1)
    
    shortening = z1/z0
    T, qs1 = geotherm_steady(depths,
                        L0*shortening,
                        shortening,
                        Ts=Ts,
                        Tlab=Tlab,
                        k=conductivity,
                        A=As,
                        hr0=hr0,)
    T1, _qs = geotherm_steady(z0/L0,
                L0*shortening,
                shortening,
                Ts=Ts,
                Tlab=Tlab,
                k=conductivity,
                A=As,
                hr0=hr0)
    plt.plot(T-273.15, depths_sc/1e3*shortening, "--", color=color,alpha=0.5,linewidth=1)

    Tts = np.zeros(100)
    zts = np.zeros(100)
    for i,t in enumerate(np.linspace(0,1,100)):
        s = 1 + (z1/z0 - 1)*t
        Tt, _q = geotherm_steady(z0/L0,
                L0*s,
                s,
                Ts=Ts,
                Tlab=Tlab,
                k=conductivity,
                A=As,
                hr0=hr0)
        Tts[i] = Tt
        zts[i] = z0*s
    
    label = setting["setting"]

    # T-z-t path
    plt.plot(Tts-273.15, zts/1.e3,'-',alpha=1,linewidth=1.2,color=color,label=label)

    # points after shortening
    plt.plot(T1-273.15, z0/1e3*shortening,'.',color=color)

    geotherm_latex_table["heading"][num_setting] = setting["setting"]
    geotherm_latex_table["As"][num_setting] = As*1.e6
    geotherm_latex_table["L0"][num_setting] = L0/1.e3
    geotherm_latex_table["z0"][num_setting] = z0/1.e3
    geotherm_latex_table["hr0"][num_setting] = hr0/1.e3
    geotherm_latex_table["qs0"][num_setting] = qs0*1.e3
    geotherm_latex_table["qs1"][num_setting] = qs1*1.e3
    geotherm_latex_table["T0"][num_setting] = T0-273.15
    geotherm_latex_table["T1"][num_setting] = T1-273.15

plt.legend()
ax1.set_ylabel("depth (km)")
ax1.set_xlabel("$T$ (°C)")

ax1.invert_yaxis()
plt.savefig(Path(output_path,"{}.{}".format("_PTt", "pdf")), metadata=pdf_metadata)
plt.savefig(Path(output_path,"{}.{}".format("_PTt", "png")))

ax1.set_ylim([0,120])
ax1.set_xlim([150,1300])
plt.savefig(Path(output_path,"{}.{}".format("_PTt_inverted", "pdf")), metadata=pdf_metadata)
plt.savefig(Path(output_path,"{}.{}".format("_PTt_inverted", "png")))

############################################
# Write LaTeX for 'tectonic setting' table #
############################################

table_body = """\\begin{{tabular}}{{cc{}}}
\\toprule
& & {} \\\\
\\midrule
$A_{{s}}$ & \\si{{\\uW\\per\\m\\cubed}} & {} \\\\
$L_{{0}}$ & \\si{{\\km}} & {} \\\\
$z_{{0}}$ & \\si{{\\km}} & {} \\\\
$h_{{r0}}$ & \\si{{\\km}} & {} \\\\
\\midrule
$q_{{s0}}$ & \\si{{\\mW\\per\\m\\squared}} & {} \\\\
$q_{{s1}}$ & \\si{{\\mW\\per\\m\\squared}} & {} \\\\
$T_{{0}}$ & \\si{{\\degreeCelsius}} & {} \\\\
$T_{{1}}$ & \\si{{\\degreeCelsius}} & {} \\\\
\\bottomrule
\\end{{tabular}}
""".format(
    "c"*len(geotherm_latex_table["heading"]),
    " & ".join(geotherm_latex_table["heading"]),
    " & ".join(["{:.2f}".format(val) for val in geotherm_latex_table["As"]]),
    " & ".join(["{:.0f}".format(val) for val in geotherm_latex_table["L0"]]),
    " & ".join(["{:.0f}".format(val) for val in geotherm_latex_table["z0"]]),
    " & ".join(["{:.1f}".format(val) for val in geotherm_latex_table["hr0"]]),
    " & ".join(["{:.0f}".format(val) for val in geotherm_latex_table["qs0"]]),
    " & ".join(["{:.0f}".format(val) for val in geotherm_latex_table["qs1"]]),
    " & ".join(["{:.0f}".format(val) for val in geotherm_latex_table["T0"]]),
    " & ".join(["{:.0f}".format(val) for val in geotherm_latex_table["T1"]]),

)

with open(Path(output_path,"_geotherms_table.tex"), "w") as fil:
    fil.writelines(table_body)

############################
# Setup initial conditions #
############################

# Create scenarios for all combinations 
# of setting, composition, and Da       
scenarios:List[InputScenario] = []
for setting in tectonic_settings:
    for da in Das:
        for comp in compositions:
            scenario = setting.copy()
            scenario["Da"] = da
            scenario["composition"] = comp
            scenarios.append(scenario)

# Get reaction/phase/endmember parameters
rxn = get_reaction(rxn_name)
phase_names, endmember_names = get_names(rxn)

I = len(rxn.phases())
_Kis = [len(rxn.phases()[i].endmembers()) for i in range(I)] # list, num EMs in each phase
K = sum(_Kis)

def reshape_C(rxn,cik:List[float])->List[List[float]]:
    c:List[List[float]] = rxn.zero_C()
    k = 0
    for i,Ki in enumerate(_Kis):
        c[i] = cik[k:k+Ki]
        k = k+Ki
    return c

def get_rho(rxn, Fi:List[float], cik:List[float], T:float, P:float)->float:
    # Calculate rho for each timestep as 1/sum_i(F_i/rho_i)
    # for which we need the endmember compositions as a vector for each phase (Cs_times)
    rhoi = rxn.rho(T,P,cik) # phase densities rho_i
    rho = 1/sum(Fi/rhoi)/10.0 # g/cm3 
    return rho

# Get equilibrated initial condition (with cache to avoid repeating work)
ic_cache = {}
def get_initial_composition(T0:float,P0:float,composition_name:str)->Tuple[List[float],List[float],float]:
    rxn = get_reaction(rxn_name)
    slug = "{:.4f}-{:.4f}-{}-{}".format(T0,P0,rxn.name(),composition_name)
    if slug in ic_cache:
        return ic_cache[slug]
    
    rxn.set_parameter("T0",Tr)

    # Get equilibrium composition from Perple_X (arbitrary [P,T])
    Fi_a, Xik_a, phii_a, cik_a = ppx_point_composition(rxn, composition_name) 
    cik_a = x2c(rxn, Xik_a) if cik_a is None else cik_a
    Fi_a = phi2F(rxn, phii_a, cik_a) if Fi_a is None else Fi_a

    # Equilibrate the reative model at initial (T0, P0) with Da=1e8
    # Set up vector of initial conditions
    u0_a = np.concatenate((Fi_a,cik_a))

    rho_a = get_rho(rxn,np.asarray(Fi_a),reshape_C(rxn,cik_a),T0,P0)*10.
    print(rho_a)
    args = (rxn,Da_eq,T0,P0,rho_a)
    sol = solve_ivp(rhs_fixed, [0, 1], u0_a, args=args, dense_output=True, method="BDF", rtol=rtol, atol=atol)
    
    Fi0 = sol.y[:I,-1] # -1 = final timestep
    cik0 = sol.y[I:I+K,-1]
    rho0 = get_rho(rxn,Fi0,reshape_C(rxn,cik0),T0,P0)*1000. # kg/m3

    ic_cache[slug] = (cik0, Fi0, rho0)
    return cik0, Fi0, rho0

# Get full initial conditions
def get_initial_state(scenario:TectonicSetting)->InputScenario:

    composition_name = scenario['composition']
    L0 = scenario["L0"]
    As = scenario["As"]
    hr0 = scenario["hr0"]
    k = scenario["k"]
    Ts = scenario["Ts"]
    Tlab = scenario["Tlab"]

    # Initial temperature
    T0, qs0 = geotherm_steady(
        z0/L0, # between 0 and 1
        L0,
        thickening=1.0, # gt 1
        Ts=Ts,
        Tlab=Tlab,
        k=k,
        A=As,
        hr0=hr0
    )

    # Final temperature
    thickening = z1/z0
    T1, qs1 = geotherm_steady(
        z0/L0, # between 0 and 1
        L0*thickening,
        thickening=thickening, # gt 1
        Ts=Ts,
        Tlab=Tlab,
        k=k,
        A=As,
        hr0=hr0
    )

    # Initial pressure
    P0 = crustal_rho * gravity * z0/1e5 # bar

    # Initial composition (equilibrium at T0,P0)
    cik0, Fi0, rho0 = get_initial_composition(T0,P0,composition_name)

    return T0, P0, cik0, Fi0, rho0, qs0, T1, qs1

# Gets all initial conditions and save to scenario
def setup_ics(scenario:TectonicSetting)->InputScenario:
    T0, P0, cik0, Fi0, rho0,qs0,T1,qs1 = get_initial_state(scenario)
    scenario["T0"] = T0
    scenario["T1"] = T1
    scenario["qs0"] = qs0
    scenario["qs1"] = qs1
    scenario["P0"] = P0
    scenario["cik0"] = cik0
    scenario["Fi0"] = Fi0
    scenario["rho0"] = rho0
    return scenario

#############################################
# RHS to equilibrate at fixed P,T           #
#############################################

def rhs_fixed(t,u,rxn,Da,T,P,rho0):

    Fi = u[:I] # phase mass fractions
    cik = u[I:I+K] # endmember mass fractions
  
    # reshape C
    C = rxn.zero_C() # object with correct shape
    Kis = 0
    for i,Ki in enumerate(_Kis):
        C[i] = cik[Kis:Kis+Ki]
        Kis = Kis+Ki

    # regularize C by taking max of C and eps
    Cs = [np.maximum(np.asarray(C[i],dtype=np.double), eps*np.ones(len(C[i]))) for i in range(len(C))]
    Cs = [np.asarray(Cs[i],dtype=np.double)/sum(Cs[i]) for i in range(len(Cs))]
    rhoi = np.array(rxn.rho(T, P, Cs)) # phase densities $\rho_i$
    V = np.sum(Fi/rhoi) # total volume
  
    # regularize F by adding eps
    Fis = np.asarray(Fi)
    Fis = Fi + eps

    # Get dimensionless Gammas from reaction
    Gammai = np.asarray(rxn.Gamma_i(T,P,Cs,Fi))
    gamma_ik = rxn.Gamma_ik(T,P,Cs,Fi)
    Gammaik = np.zeros(K)
    sKi = 0
    for i in range(I):
        for k in range(_Kis[i]):
            Gammaik[sKi+k] = gamma_ik[i][k]
        sKi += _Kis[i]
    
    # Calculate reactive mass flux (scaled Gammas)
    du = np.zeros(u.shape)
    sKi = 0
    for i in range(I):
        du[i] = Da*rho0*Gammai[i]*V
        for k in range(_Kis[i]):
            GikcGi = Gammaik[sKi+k] - C[i][k]*Gammai[i]
            du[I+sKi+k] = Da*rho0*GikcGi*V/Fis[i]
        sKi += _Kis[i]
    return du

#############################################
# Define RHS of differential system of eqns #
#############################################

def rhs(t,u,rxn,scale,thermal):
    
    # Extract variables
    Fi = u[:I] # phase mass fractions
    cik = u[I:I+K] # endmember mass fractions

    h0 = scale["h"]
    h0 = scale["h"]
    rho0 = scale["rho"]
    Da = scale["Da"]

    L0 = thermal["L0"]
    z0 = thermal["z0"]
    As = thermal["As"]
    hr0 = thermal["hr0"]
    k = thermal["k"]
    Ts = thermal["Ts"]
    Tlab = thermal["Tlab"]

    # limiting depth to some value (e.g. 200 km) required here. If python tries to take a very large timestep
    # we will be out of bounds for the thermodynamic database
    z_t = min(200e3, z0 + t*h0)

    shortening_t = z_t/z0
    P_t = crustal_rho * gravity * z_t / 1e5 # bar
    T_t, q_s = geotherm_steady(z0/L0,
        L0*shortening_t,
        shortening_t,
        Ts=Ts,
        Tlab=Tlab,
        k=k,
        A=As,
        hr0=hr0) # K
    
    # reshape C
    C = rxn.zero_C() # object with correct shape
    Kis = 0
    for i,Ki in enumerate(_Kis):
        C[i] = cik[Kis:Kis+Ki]
        Kis = Kis+Ki

    # regularize C by taking max of C and eps
    Cs = [np.maximum(np.asarray(C[i]), eps*np.ones(len(C[i]))) for i in range(len(C))]
    Cs = [np.asarray(Cs[i])/sum(Cs[i]) for i in range(len(Cs))]

    rhoi = np.array(rxn.rho(T_t, P_t, Cs)) # phase densities $\rho_i$
    V = np.sum(Fi/rhoi) # total volume
    
    # regularize F by adding eps
    Fis = np.asarray(Fi) + eps
    Fis = Fi + eps

    # Get dimensionless Gammas from reaction
    Gammai = np.asarray(rxn.Gamma_i(T_t,P_t,Cs,Fi))
    gamma_ik = rxn.Gamma_ik(T_t,P_t,Cs,Fi)
    Gammaik = np.zeros(K)
    sKi = 0
    for i in range(I):
        for k in range(_Kis[i]):
            Gammaik[sKi+k] = gamma_ik[i][k]
        sKi += _Kis[i]
    
    # Calculate reactive mass flux (scaled Gammas)
    du = np.zeros(u.shape)
    sKi = 0
    for i in range(I):
        du[i] = Da*rho0*Gammai[i]*V
        for k in range(_Kis[i]):
            GikcGi = Gammaik[sKi+k] - C[i][k]*Gammai[i]
            du[I+sKi+k] = Da*rho0*GikcGi*V/Fis[i]
        sKi += _Kis[i]

    return du

############################################
# Function that runs scenarios in parallel #
############################################

def run_experiment(scenario:InputScenario)->OutputScenario:
    Da = scenario["Da"]
    L0 = scenario["L0"]
    As = scenario["As"]
    hr0 = scenario["hr0"]
    cik0 = scenario["cik0"]
    Fi0 = scenario["Fi0"]
    rho0 = scenario["rho0"]
    k = scenario["k"]
    Ts = scenario["Ts"]
    Tlab = scenario["Tlab"]

    rxn = get_reaction(rxn_name)

    # Set reaction's characteristic Arrhenius temperature (T_r)
    rxn.set_parameter("T0",Tr)

    # Set up vector of initial conditions
    u0 = np.empty(I+K) # [...I phases, ...K endmembers]
    u0[:I] = Fi0 # intial phase mass fractions
    u0[I:I+K] = cik0 # initial endmember mass fractions

    scale = {"rho":rho0, "h":h0, "Da":Da}
    thermal = {"L0":L0, "z0":z0, "As":As,"hr0":hr0,"k":k,"Ts":Ts,"Tlab":Tlab}
    args = (rxn, scale, thermal)

    # Solve IVP using BDF method
    sol = solve_ivp(rhs, [0, end_t], u0, args=args, dense_output=True, method="BDF", rtol=rtol, atol=atol, events=None)
    
    # resample solution
    times = np.linspace(0,end_t,1000)
    y = sol.sol(times)

    Fi_times  = y[:I].T # vector for each timestep
    cik_times = y[I:I+K].T # 2d ragged array for each timestep
    Cs_times = [reshape_C(rxn,cik) for cik in cik_times] # vector for each timestep

    # Back-calculate depth, T, P, and rho
    depth_m_times = z0 + times*h0
    shortening_times = depth_m_times/z0
    T_times = [geotherm_steady(
        z0/L0,
        L0*shortening,
        shortening,
        Ts=Ts,
        Tlab=Tlab,
        k=conductivity,
        A=As,
        hr0=hr0)[0] for shortening in shortening_times]
    P_times = crustal_rho*gravity*depth_m_times/1e5 # bar
    
    # Calculate rho for each timestep as 1/sum_i(F_i/rho_i)
    # for which we need the endmember compositions as a vector for each phase (Cs_times)
    rho_times = [1/sum(Fi_times[idx]/rxn.rho(T_times[idx], P_times[idx], Cs_times[idx]))/10. for idx,_ in enumerate(times)]
    print("{} P_end = {:.2f} Gpa. T_end = {:.2f} K. DA = {}. Used {:n} steps.".format(sol.message,P_times[-1]/1e4,T_times[-1],Da,len(sol.t)))

    scenario["T"] = np.asarray(T_times) # K
    scenario["P"] = np.asarray(P_times) # bar
    scenario["rho"] = np.asarray(rho_times) # g/cm3
    scenario["Fi"] = Fi_times # phase mass fractions
    scenario["cik"] = cik_times # endmember mass fractions
    scenario["Xik"] = np.asarray([rxn.C_to_X(c) for c in Cs_times], dtype="object") # endmember mol. fractions
    scenario["z"] = depth_m_times
    scenario["time"] = times # 
    return scenario

#####################################
# Deal with result loading & saving #
#####################################
scenarios_out = None
if load_output:
    print("Looking for pickle file {}".format(pickle_path))
    try:
        with open(pickle_path, 'rb') as pickle_file:
            scenarios_out = pickle.load(pickle_file)
            print("Successfully loaded output")
    except:
      sys.stdout.write("\n"+bcolors.WARNING+"WARNING: Unable to load file {}".format(pickle_path)+bcolors.ENDC)
      sys.stdout.write("\n"+bcolors.WARNING+"Re-calculate all model scenarios (may take some time)"+bcolors.ENDC)
      sys.stdout.write("\n"+bcolors.OKCYAN+"Continue? [Y/n]: "+bcolors.ENDC)
      yes = {'yes','y', 'ye', ''}
      no = {'no','n'}

      choice = input().lower()
      if choice not in yes:
        print("Quitting...")
        quit()
      load_output = False
      save_output = True

if scenarios_out is None:
    print("Preparing to run {} scenarios".format(len(scenarios)))
    scenarios_in = [setup_ics(s) for s in scenarios]

    # run for varying damkhoeler numbers
    with Pool(num_processes) as pool:
        # blocks until all finished
        scenarios_out = pool.map(run_experiment, scenarios_in)

if save_output:
    with open(pickle_path, 'wb') as pickle_file:
        pickle.dump(scenarios_out, pickle_file)

###################
# Post processing #
###################

# this needs to be in global scope
bin_widths = {"1Myr":1*Myr, "100kyr":100*kyr, "50kyr":50*kyr, "10kyr":10*kyr}

for out in scenarios_out:
    # Process each scenario output
    rho = np.array(out["rho"]) # g/cm3
    depth_m = out["z"]
    depth_lab = depth_m/z0 * out["L0"]
    T = out["T"] # K
    t = out["time"] # unitless
    P = out["P"] # bar
    rho_pyrolite = model_pyrolite_rho_gcc(T,P)
    max_rho = np.nanmax(rho)
    max_rho_py = np.nanmax(rho_pyrolite)

    # find the critical depth (pressure, temperature, etc)
    if max_rho > max_rho_py:
        # find the greatest depth at which rho exceeds pyrolite
        critical_indices = [i for i,r in enumerate(rho) if r > rho_pyrolite[i]]
        if len(critical_indices) == 0:
            # never goes critical
            critical_depth = np.nan # 108.e3 # approx. coesite transition
            critical_pressure = np.nan # 30.e3 # bar
            critical_temperature = np.nan # T[-1] # K
            critical_time = np.nan # end_t+1
        elif len(critical_indices) == len(rho):
            # always critical
            critical_depth = depth_m[0]
            critical_pressure = P[0] # bar
            critical_temperature = T[0] # K
            critical_time = t[0]
        else:
            # interesting 
            first_critical_index = critical_indices[0]
            critical_index = first_critical_index
            critical_depth = depth_m[critical_index]
            critical_pressure = P[critical_index] # bar
            critical_temperature = T[critical_index] # K
            critical_time = t[critical_index]
    else:
        critical_depth = np.nan # 108.e3 # approx. coesite transition
        critical_pressure = np.nan # 30.e3 # bar
        critical_temperature = np.nan # T[-1] # K
        critical_time = np.nan # end_t+1
    out["critical_depth"] = critical_depth
    out["critical_pressure"] = critical_pressure
    out["critical_temperature"] = critical_temperature
    out["critical_time"] = critical_time

    # Average density contrast of any unstable root
    delta_rho = rho - rho_pyrolite
    delta_rho[delta_rho < 0] = np.nan
    effective_delta_rho = np.nanmean(delta_rho)
    out["effective_delta_rho"] = effective_delta_rho
    out["max_rho"] = max_rho
    out["max_delta_rho"] = max_rho - max_rho_py if max_rho - max_rho_py > 0 else float("NaN")

    # Densification rate
    time = np.linspace(0,end_t, rho.size) * t0 # seconds
    time_Myr = time/Myr
    densification_rate = np.diff(rho*1000)/np.diff(time_Myr) # kg/m3/Myr
    densification_rate = np.insert(densification_rate, 0, 0.)
    out["densification_rate"] = densification_rate
    out["time_Myr"] = time_Myr
    for i,phase in enumerate(phase_names):
        phase_mis = out["Fi"][:,i] # wt%
        if(phase=='Feldspar'):
            plag_frac = phase_mis.copy()
            # find the index at which plagioclase drops below 1 wt%
            plag_out_indices = [i for i,X in enumerate(plag_frac) if X < 0.025]
            if len(plag_out_indices) == 0:
                out['plag_out_depth'] = np.nan 
                out['plag_out_pressure'] = np.nan # bar
                out['plag_out_temperature'] = np.nan # K
            else:
                first_plag_out = plag_out_indices[0]
                plag_out_index = first_plag_out
                out['plag_out_depth'] = depth_m[plag_out_index]
                out['plag_out_pressure'] = P[plag_out_index] # bar
                out['plag_out_temperature'] = T[plag_out_index] # K

    # Max densification rate, binned and averaged over time bins,
    # limited to the P-T space of the plag-out reaction (not the garnet-in reaction)
    for bin_width_string, bin_width in bin_widths.items():
        bins = np.arange(0., t0, int(bin_width))
        digitized_t = np.digitize(time, bins) # assigns the index of a bin to each point
        plag_out_mask = P/1.e4 < 0.5 + 1./1000.*(T-273.15-300.) # NOTE: 'True' signifies masking OUT (nan)
        dens_rate_masked = ma.masked_array(densification_rate,mask=plag_out_mask)
        bin_means = [dens_rate_masked[digitized_t == i].mean() for i in range(len(bins))]
        out["max_densification_rate_"+bin_width_string] = np.nanmax(bin_means)

#############################
# Write summary data to CSV #
#############################

with open(Path(output_path,'_summary.csv'),'w') as csvfile:
    fieldnames = [
        'setting',
        'L0',
        'z0',
        'z1',
        'As',
        'hr0',
        'k',
        'Ts',
        'Tlab',
        'Da',
        'composition',
        'T0',
        'T1',
        'P0',
        'rho0',
        'critical_depth',
        'critical_pressure',
        'critical_temperature',
        'plag_out_pressure',
        'plag_out_temperature',
        'plag_out_depth',
        'qs0',
        'qs1',
        'effective_delta_rho',
        'max_delta_rho',
        'max_rho'
    ]
    for key,val in bin_widths.items():
        fieldnames.append('max_densification_rate_'+key)
    writer = csv.DictWriter(csvfile, fieldnames,extrasaction='ignore')
    writer.writeheader()
    for out in scenarios_out:
        writer.writerow(out)
    print("wrote CSV")

#############################
# Rayleigh--Taylor analysis #
#############################

selected_compositions = ["sammon_2021_lower_crust","hacker_2015_md_xenolith","mackwell_1998_maryland_diabase"]
fluid_weakening = [1, 0.5, 0.25, 0.1] # Weaken B to account for fluids (base eclogite rheology is dry)
for f in fluid_weakening:
    for _da in [10,300,1000,3000,10000]:
        # Setup figure for Rayleigh-Taylor analysis by composition
        fig1 = plt.figure(figsize=(3.75, 3.5))
        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
        # y axis is time (yr) in log scale
        # x axis is layer thickness (km), linear
        plt.gca().set_ylim([1e5,5e7])
        plt.gca().set_yscale("log")
        plt.gca().set_xlim([0,40])

        v0_yr = v0 * yr # converted from m/s to m/yr
        t_growth_log = np.logspace(5,8,num=200) # yr
        h_growth_log = t_growth_log*v0_yr/1e3 # km

        plt.plot(h_growth_log, t_growth_log,'k-')
        plt.plot(h_growth_log*0.7, t_growth_log,'k--',alpha=0.7,linewidth=0.5)
        plt.plot(h_growth_log*1.5,t_growth_log,'k--',alpha=0.7,linewidth=0.5)
        plt.gca().set_prop_cycle(plt.cycler("linestyle", ['-','--',':']))

        outs_c = sorted([o for o in scenarios_out 
            if (o["composition"] in selected_compositions and o["Da"]==_da)
        ], key=lambda o: (o["composition"],o["setting"]))

        for i, obj in enumerate(outs_c):
            # get values from simulation
            composition = obj["composition"]
            color = color_by_composition.get(composition, "black")
            setting = obj["setting"]
            T = obj["T"]
            P = obj["P"]
            max_temp = obj["T1"]
            Da = obj["Da"]
            critical_h = obj["critical_depth"]
            rho_pyrolite = model_pyrolite_rho_gcc(T,P)*1000. # k/m3
            densities = np.array(obj["rho"]) * 1000.

            # Define "root" as negatively buoyant portion (drho > 0)
            drho = densities - rho_pyrolite
            is_root = drho > 0.
            max_drho = max(drho)
            if(max_drho <= 0):
                continue
            
            # pull out drho and T of just the root
            root_drho = drho[is_root]
            root_T = T[is_root]
            
            # set up thickness variable
            # note: because we are shifting portions of arrays around,
            # both the length of these arrays and the size of dt have to match how the solution 'y' was calculated
            t_growth_yr = obj["time"]*t0/yr # 0 to 50e6 years
            h = t_growth_yr/1e3  # 0 to 50e3 meteres

            # allow root to thicken to 40 km, extending as needed with the final value from the simulation
            root_drho_extended = np.ones(h.size)*root_drho[-1]
            root_drho_extended[:root_drho.size] = root_drho
            T_extended = np.ones(h.size)*root_T[-1]
            T_extended[:root_T.size] = root_T

            # Eclogite (Jin et al. 2001), following Molnar & Garzione, Zieman
            Rgas = 8.3145
            g = 9.81
            A_Mpa = 10.**3.3
            Q = 480.e3
            n = 3.4
            
            # Wet olivine
            # A_Mpa = 1.9e3
            # Q = 420e3
            # n = 3.

            A = A_Mpa*(1e6)**(-n) # Pa^-3.4 s^-1
            F = 3.**(-(n+1.)/2./n)*(2)**(1./n) # convert imposed strain fields in lab to a general geometry
            B = f*F*(A)**(-1./n)*np.exp(Q/(n*Rgas*T_extended)) # Pa s = kg/m/s

            # High-stress version of B (Evans & Goetze 1979, Molnar & Jones 2004)
            eta0 = 5.7e11 # 1/s
            sigma0 = 8.5e9 # Pa
            Ha = 525.e3 # J/mol/K
            E = 1.e-14 # sqrt of 2nd invariant strain rate (1/s)
            B_Peierls = E**((n-1)/n) * sigma0 / (E*np.sqrt(3.0)) * (1. - np.sqrt((Rgas*T_extended)/(Ha) * np.log((np.sqrt(3.0)*eta0)/(2.0*E))) )
            
            # Effective viscosity
            B_eff = B
            B_eff[T_extended<1000] = np.minimum(B[T_extended<1000], B_Peierls[T_extended<1000])

            # Time-dependent, moving average of drho and B_eff as root grows
            avg_drho = np.array([np.average(root_drho_extended[:i+1]) for i, r in enumerate(root_drho_extended)])
            avg_B_eff = np.array([np.average(B_eff[:i+1]) for i, b in enumerate(B_eff)])

            # Growth rate and perturbation factors, Jull & Kelemen (2001)
            Cp = 0.66 # strong layer, L>>h, following Zieman et al.
            Zp0 = 0.33 # Assume 33% for initial perturbation amplitude, follows Zieman et al.

            Timescale = (avg_B_eff/(2.*avg_drho*g*h))**n # Eq 7, timescale in seconds
            tbp0 = ((n/Cp)**n)*((Zp0)**(1-n))/(n-1) # Eq 12, dimensionless time for 100% deflection

            exx = 1e-14 # horizontal strain rate, Behn et al. 2007, Zieman et al. 2023
            epxx = exx * Timescale # by Eq. 8, dimensionless, approx 10
            epxx0 = 1e-18 * Timescale # by Eq. 8, dimensionless, approx 1e-3
            dtpdep = -0.5 # log units

            exponent = np.double(-epxx/tbp0*dtpdep)
            tbp = tbp0*(epxx/epxx0)**exponent # instability time, dimensionless

            tb_yr = tbp*Timescale/yr # instability time, years

            plt.figure(fig1)
            plt.plot(h[h>z1-critical_h]/1e3, tb_yr[h>z1-critical_h], '-', linewidth=(T_extended[-1]/1273.15)**2, color=color,alpha=0.25)
            plt.plot(h[h<=z1-critical_h]/1e3, tb_yr[h<=z1-critical_h], linewidth=(T_extended[-1]/1273.15)**2, color=color,label=composition)
            plt.plot(z1-critical_h, (z1-critical_h)/v0/yr, 'o',color=color)  

            intersection_idx = np.argwhere(np.diff(np.sign(tb_yr - t_growth_yr))).flatten()
            intersection_idx = intersection_idx[intersection_idx > 0]
            if(intersection_idx.size):
                intersection_idx = intersection_idx[0]
                final_h = h[intersection_idx]
                final_tb = tb_yr[intersection_idx]
                final_T = T_extended[intersection_idx]
            else:
                final_h = np.nan
                final_tb = np.nan
                final_T = np.nan

            plt.plot(final_h/1.e3, final_tb, 'o', color=color,markersize=4)

            if(_da==10000):
                obj["final_h_{}".format(f)] = final_h
                obj["final_depth_{}".format(f)] = final_h + obj["critical_depth"]
                obj["final_T_{}".format(f)] = final_T
                obj["final_t_{}".format(f)] = final_tb + obj["critical_time"]
        
        plt.figure(fig1)
        plt.savefig(Path(output_path,"_instability.Da{}.f{}.{}".format(_da,f,"pdf")), metadata=pdf_metadata)
        plt.savefig(Path(output_path,"_instability.Da{}.f{}.{}".format(_da,f,"png")))
        plt.close(fig1)

###########################
# Plot max. stable depths #
###########################

fig = plt.figure(figsize=(20,7.5))
axes = fig.subplot_mosaic([selected_compositions])

# Invert y axis because it represents depth
[ax.invert_yaxis() for label,ax in axes.items()]
[ax.set_ylim([z1/1.e3,z0/1.e3]) for label,ax in axes.items()]
[ax.tick_params(width=0.4) for label,ax in axes.items()]
for axis in ['top','bottom','left','right']:
    [ax.spines[axis].set_linewidth(0.25) for label,ax in axes.items()]

[ax.set_xlabel("Temperature (°C)")for label,ax in axes.items()]
[ax.set_ylabel("Critical depth (km)")for label,ax in axes.items()]

for comp in selected_compositions:
    ax = axes[comp]
    color = color_by_composition.get(comp, "black")

    for _da in Das:
        outs_c_da = sorted([out for out in scenarios_out if (out["composition"] == comp) and out["Da"]==_da], key=lambda out: out["setting"],reverse=True)
        crit_T = np.array([o["critical_temperature"] for o in outs_c_da])
        crit_z = np.array([o["critical_depth"] for o in outs_c_da])
        s = ax.plot(crit_T-273.15, crit_z/1.e3,color=color,linewidth=0.75)
        if _da==10000:
            for f in fluid_weakening:
                final_z = np.array([o["final_depth_{}".format(f)] for o in outs_c_da])
                final_T = np.array([o["final_T_{}".format(f)] for o in outs_c_da])
                ax.plot(final_T-273.15, final_z/1.e3, 'r--',alpha=np.sqrt(f))
        for obj in outs_c_da:
            ax.plot(obj["T"]-273.15, obj["z"]/1.e3, color='#888888',linewidth=0.25, label=obj['setting'])

plt.savefig(Path(output_path,"{}.{}".format("_critical", "pdf")), metadata=pdf_metadata)
plt.savefig(Path(output_path,"{}.{}".format("_critical", "png")))
plt.close(fig)

##########################################################
# Summary plots for every composition & tectonic setting #
##########################################################
    
for composition in compositions:
    for tectonic_setting in tectonic_settings:
        setting = tectonic_setting["setting"]

        # Find outputs with same setting & same composition,
        # and sort by Da
        outs_c = sorted([out for out in scenarios_out if (out["composition"] == composition and out["setting"] == setting)], key=lambda out: out["Da"])
        
        # grab 'constants' from first output
        base = outs_c[0] 
        T = base["T"]
        P = base["P"]
        rho0 = base["rho0"]
        L0 = base["L0"]

        _Das = [o["Da"] for o in outs_c]
        
        # Back-calculate reaction coefficient from Da
        r0 = [da*rho0/t0 for da in _Das] # Gamma0 (kg/m3/s)
        S0 = 6000 # 1/m
        reaction_rate_per_surface = [r/S0 for r in r0] # r0 (kg/m2/s)
        reaction_rate_per_surface_gcm = [r*1000/100/100*yr for r in reaction_rate_per_surface] # g/cm2/yr

        # Setup subplots and figure
        num_subplots = 3 + len(phase_names) + 2
        subplot_mosaic = [["rho","rho"] + phase_names + ["An", "Jd", "T"]] # axes are named
        fig = plt.figure(figsize=(3*num_subplots,12))
        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
        axes = fig.subplot_mosaic(subplot_mosaic)

        # Invert y axis because it represents depth
        [ax.invert_yaxis() for label,ax in axes.items()]

        # Setup color cycling through all Das
        num_lines = len(_Das)
        greys = plt.cm.get_cmap("Greys")
        axes["rho"].set_prop_cycle(plt.cycler("color", greys(np.linspace(0.2, 1, num_lines))))
        axes["An"].set_prop_cycle(plt.cycler("color", greys(np.linspace(0.2, 1, num_lines))))
        axes["Jd"].set_prop_cycle(plt.cycler("color", greys(np.linspace(0.2, 1, num_lines))))

        # Plot density vs depth
        for i, obj in enumerate(outs_c):
            ax = axes["rho"]
            ax.plot(obj["rho"], obj["z"])

        rho_pyrolite = model_pyrolite_rho_gcc(T,P)
        ax.plot(rho_pyrolite/10, outs_c[0]["z"], "r:")
        ax.legend(["$Da = ${:.1e}".format(d) for d in _Das] +  (["pyrolite"]), loc="upper right")

        ax.set_ylabel("Depth (km)")
        ax.set_xlabel("Density")
        ax.set_xlim([2.8, 3.8])

        # Plot temperature vs depth
        ax = axes["T"]
        ax.plot(T-273.15, base["z"]/1e3, linewidth=2)
        ax.set_xlabel("T (°C)")
        ax2 = ax.twinx()
        ax2.plot(T-273.15, P/1e4, alpha=0)
        ax2.invert_yaxis()
        ax2.set_ylabel("P (GPa)")

        # Plot all phase mis (mass fractions)
        cmaps = [plt.cm.get_cmap(name) for name in ["Blues", "YlOrBr", "Greens","Reds","Purples","copper_r"]]
        for i,phase in enumerate(phase_names):
            ax = axes[phase]
            cmap = cmaps[i]
            ax.set_prop_cycle(plt.cycler("color", cmap(np.linspace(0.2, 1, num_lines))))
            ax.set_xlim([0., z1/1.e3])
            ax.set_ylabel(None)
            ax.set_xlabel("{} (wt%)".format(phase))
            ax.set_xticks(np.arange(0,110,10))
            ax.set_xticklabels([None, None, 20, None, 40, None, 60, None, 80, None, None])
            ax.set_yticklabels([])
            for j, obj in enumerate(outs_c):
                ax.plot(obj["Fi"][:,i]*100, obj["z"])

        # get Plagioclase anorthite composition
        ax = axes["An"]
        ax.set_xlim([0., 1.0])
        ax.set_ylabel(None)
        ax.set_yticklabels([])
        ax.set_xlabel("$X_{\mathrm{An}}$")

        # get Cpx jadeite content
        ax = axes["Jd"]
        ax.set_xlim([0., 1.0])
        ax.set_ylabel(None)
        ax.set_yticklabels([])
        ax.set_xlabel("$X_{\mathrm{Jd}}$")

        # Plot X_An and X_Jd vs depth
        for j, obj in enumerate(outs_c):
            XAn = [x[3][0] for x in obj["Xik"]]
            XJd = [x[0][4] for x in obj["Xik"]]

            axes["An"].plot(XAn, obj["z"])
            axes["Jd"].plot(XJd, obj["z"])

        fig.suptitle("{}, {}, $v_0=${:.1f} km/Myr, $S0=${} 1/m".format(composition.capitalize().replace("_"," "), setting.replace("_",", "), v0/1e3*yr*1e6, S0),y=0.9)
        plt.savefig(Path(output_path,"{}.{}.{}".format(setting,composition,"pdf")), metadata=pdf_metadata)
        plt.savefig(Path(output_path,"{}.{}.{}".format(setting,composition,"png")))
        plt.close(fig)


#####################################################
# For each setting, plot each composition's density #
#####################################################
selected_compositions = [
    "sammon_2021_deep_crust",
    "sammon_2021_lower_crust",
    "hacker_2015_md_xenolith",
    "mackwell_1998_maryland_diabase"
]        
for tectonic_setting in tectonic_settings:

        setting = tectonic_setting["setting"]
        outs_c = sorted([out for out in scenarios_out if (out["composition"] in selected_compositions and out["setting"] == setting)], key=lambda out: out["Da"])
        
        # grab 'constants' from first output
        base = outs_c[0] 
        T = base["T"]
        P = base["P"]
        # Setup figure for rho-profile mosaic
        fig = plt.figure(figsize=(6.5,3.0))
        axes = fig.subplot_mosaic([selected_compositions])
        plt.tight_layout(pad=0.075, w_pad=0.075, h_pad=.075)
        [ax.set_prop_cycle(plt.cycler("color", greys(np.linspace(0.2, 1, num_lines)))) for label,ax in axes.items()]

        # Invert y axis because it represents depth
        [ax.invert_yaxis() for label,ax in axes.items()]
        [ax.set_xlim([2.8,3.5]) for label,ax in axes.items()]
        [ax.set_xticks([2.8,3.0,3.2,3.4]) for label,ax in axes.items()]
        [ax.set_ylim([z1/1.e3,z0/1.e3]) for label,ax in axes.items()]
        [ax.tick_params(width=0.4) for label,ax in axes.items()]
        for axis in ['top','bottom','left','right']:
            [ax.spines[axis].set_linewidth(0.25) for label,ax in axes.items()]
  
        rho_pyrolite = model_pyrolite_rho_gcc(T,P)

        for i, obj in enumerate(outs_c):
            ax = axes[obj["composition"]]
            color=color_by_composition.get(obj["composition"],"black")
            linewidth = 0.5 if obj["Da"] == np.max(Das) else 0.25
            alpha = 1 if obj["Da"] == np.max(Das) else 0.8
            ax.plot(obj["rho"], obj["z"]/1e3, color=color,linewidth=linewidth)

            if obj["composition"] != selected_compositions[0]:
                ax.set_yticks([])
                ax.set_yticklabels([])
            
            if obj["Da"] == 1:
                ax.plot(rho_pyrolite, obj["z"]/1e3, "k", dashes=[10,4], linewidth=0.35)

        plt.savefig(Path(output_path,"_collage.{}.{}".format(setting,"pdf")), metadata=pdf_metadata)
        plt.savefig(Path(output_path,"_collage.{}.{}".format(setting,"png")))
        plt.close(fig)

