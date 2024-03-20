import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir, 'tcg_slb','python'))

from python.tcg import x2c,phi2F,get_reaction,composition_to_label, custom_solve
from python.perplex import ppx_point_composition, ppx_rho_interpolator
import numpy as np
from matplotlib import pyplot as plt
from pathlib import Path
from tcg_slb.phasediagram.scipy import ScipyPDReactiveODE
from tcg_slb.phasediagram.base import GPa2Bar

from multiprocessing import Pool
import multiprocessing as mp
from scipy.interpolate import griddata

SMALL_SIZE = 9
MEDIUM_SIZE = 11
BIGGER_SIZE = 13

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE, titleweight="medium")     # fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=MEDIUM_SIZE)  # fontsize of the figure title
#plt.rc('axes', titley=1.0)    # y is in axes-relative coordinates.
#plt.rc('axes', titlepad=-16)  # pad is in points...

### ======================= INPUTS ============================
reference = "parallel_pd"
composition = "hacker_2015_md_xenolith"
rxn_name = "eclogitization_2024_stx21_rx"

# number of nodes in each dimension (T,P) = (x,y)
nP = 80
nT = 80

# end time of reactions
# override with -e arg
end_t = 1

# limit the maximum number of steps taken by the solver
# defaults to infinity, which can exhaust memory if solver doesn't converge
# a default of 4e3 seems to work well
max_steps = 4e3

# regularization parameter for compositions
# sets the minimum value for cik
eps = 1.e-5 # default, 1.e-2

# TODO: these numbers seem to work very well with eps = 1e-5??
rtol = 1.e-5 # relative tolerance, default 1e-5
atol = 1.e-9 # absolute tolerance, default 1e-9

# Damkhoeler number - override with -d arg
Da = 1e6
Tr = 5500.+273.15 # reaction's characteristic temperature (T_r)

Pmin, Pmax = [0.5, 2.5] # Gpa
Tmin, Tmax = [300+273.15, 1300+273.15] # K

# Account for dense oxides not included in SLB database
oxide_density_gcc = 0.03

# override with -n arg
num_processes = mp.cpu_count()-1

# Plotting options
phasetol = 1.e-5 # phases less than this cutoff don't plot, default 1.e-2
density_levels = np.arange(2600, 3820, 20)/1000.
density_ticks = np.asarray([2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8])

T_ticks = np.arange(Tmin-273.15, Tmax-273.15+200, 200)
T_tick_labels = [""] + ["{:d}".format(int(T)) for T in T_ticks[1:-1]] + [""]
T_limits = [Tmin-273.15, Tmax-273.15]

P_ticks = np.arange(Pmin, Pmax+0.5, 0.5)

highlight_densities = np.asarray([3.33] )# /10
highlight_diffs= np.asarray([0.] )

density_cmap = plt.get_cmap("Spectral_r")
diff_cmap = plt.get_cmap("bwr_r")
stime_cmap = plt.get_cmap("bwr")
stime_cmap.set_bad(color='black') # make NaNs show up as black
variance_cmap = plt.get_cmap("Blues_r")
contour_kwargs = {
    "colors": "#333333",
    "linewidths": 1,
    "linestyles": "solid",
    "linewidths": 1,
    "alpha": 0.6
}
contour_highlight_kwargs = {
    "alpha":0.6,
    "colors":"black",
    "linewidths":1,
}
imshow_kwargs = {
    "extent": (Tmin-273.15,Tmax-273.15,Pmin,Pmax),
    "origin":"lower",
    "aspect":"auto",
}

# This custom formatter removes trailing zeros, e.g. "1.0" becomes "1", and
# then adds a percent sign.
def diff_fmt(x):
    s = f"{x:.1f}"
    if s.endswith("0"):
        s = f"{x:.0f}"
    return rf"{s}"

def rel_err_fmt(x):
    s = f"{x:.1f}"
    return rf"{s}"

fig_filetypes = ["png", "pdf"]

# ===================================================


# ============= Parse arguments for CLI =============

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--composition")
    parser.add_argument("-e", "--end_t")
    parser.add_argument("-r", "--rxn_name")
    parser.add_argument("-n", "--num_processes")
    parser.add_argument("-d", "--damk")
    args = parser.parse_args()

    if args.composition is not None:
        print("Using composition {}".format(args.composition))
        composition = args.composition
    if args.end_t is not None:
        print("Using end time {}".format(args.end_t))
        end_t = float(args.end_t)
    if args.rxn_name is not None:
        print("Using reaction {}".format(args.rxn_name))
        rxn_name = args.rxn_name
    if args.num_processes is not None:
        num_processes = int(args.num_processes)
    if args.damk is not None:
        Da = float(args.damk)

#====================================================

# get reaction object
rxn = get_reaction(rxn_name)
rxn.set_parameter("T0", Tr)
Fi0, Xik0, phii0, cik0 = ppx_point_composition(rxn, composition)

cik0 = x2c(rxn, Xik0) if cik0 is None else cik0
Fi0 = phi2F(rxn, phii0, cik0, eps) if Fi0 is None else Fi0

T_range = np.linspace(Tmin, Tmax, nT)
P_range = np.linspace(Pmin, Pmax, nP)
P_g, T_g  = np.meshgrid(P_range, T_range, indexing="ij")

rho_g = np.zeros(P_g.shape)
rho_g_oxides = np.zeros(P_g.shape)
stime_g = np.zeros(P_g.shape)
variance_g = np.zeros(P_g.shape)
phases_g = [["" for j in range(nT)] for i in range(nP)]

# function to run in parallel
# solves an ODE for each point in the grid given i, j indices
def task(args):
    i, j = args
    P = P_range[i]
    T = T_range[j]

    ode = ScipyPDReactiveODE(rxn)

    custom_solve(ode,T,GPa2Bar(P),Fi0,cik0,end_t,Da=Da,eps=eps,rtol=rtol,atol=atol)
    odephasenames, phaseabbrev = ode.final_phases(phasetol)
    phases = "+".join(phaseabbrev) if ode.sol.status == 0 else ""
    rho = ode.final_rho() if ode.sol.status == 0 else np.nan
    stime = ode.stime
    variance = odephasenames.size
    return rho/10., phases, stime, variance, i, j

# Solve individual points in parallel
with Pool(num_processes, maxtasksperchild=12) as pool:
    arglist = [[i,j] for i,P in enumerate(P_range) for j,T in enumerate(T_range)]
    sols = pool.map(task, arglist)

# Collect individual points back into grids
for rho, phases, stime, variance, i , j in sols:
    rho_g[i][j] = rho # g/cm3
    rho_g_oxides[i][j] = rho + oxide_density_gcc
    phases_g[i][j] = phases
    stime_g[i][j] = stime
    variance_g[i][j] = variance


# Find nans in density (solution timed out)
x,y = np.indices(rho_g.shape)
rhonan = np.isnan(rho_g)

# For points that timed out, give them a solution time of NaN
stime_g[rhonan] = np.nan

# Interpolate NaNs in density
rho_g[rhonan] = griddata((x[~rhonan], y[~rhonan]),rho_g[~rhonan],(x[rhonan], y[rhonan]))
rho_g_oxides[rhonan] = griddata((x[~rhonan], y[~rhonan]),rho_g[~rhonan],(x[rhonan], y[rhonan]))

# Generate unique phase assemblage indices
uniquestrs = sorted(list(set([phstr for phstrl in phases_g for phstr in phstrl if phstr != ""])))
def index(ls,v):
    try:
        i = ls.index(v)
    except ValueError:
        i = -1
    return i
phaseis_g = np.asarray([[index(uniquestrs,phstr) for phstr in phasestrr]  for phasestrr in phases_g])

# Get pyrolite density grid
interp_py = ppx_rho_interpolator("xu_2008_pyrolite", "1000kg/m3")
rho_pyrolite_g = interp_py((T_g, P_g))

# Create folders if needed
outputPath = Path("output",reference,composition,rxn_name)
outputPath.mkdir(parents=True, exist_ok=True)

def save_current_fig_as(name):
    _name = name[1:] if name[0]=="." else name
    for ext in fig_filetypes:
        _ext = ext[1:] if ext[0]=="." else ext
        plt.savefig(Path(outputPath,"{}.{}".format(_name, _ext)))

def plot_phase_labels(ax):
    fxgrid = T_g.flatten()-273.15
    fygrid = P_g.flatten()
    fphaseis = phaseis_g.flatten()

    for i,ph in enumerate(uniquestrs):
        if i in fphaseis:
            indices = fphaseis==i
            ax.text(np.mean(fxgrid[indices]), np.mean(fygrid[indices]),ph,fontsize=12,ha='center',va='center',weight="bold")

# Plot density

fig = plt.figure(figsize=(12,14))
axi = fig.add_subplot(1,1,1)

s = axi.contourf(T_g-273.15,P_g,rho_g_oxides,levels=density_levels,alpha=0.75,cmap=density_cmap)
axi.contour(T_g-273.15,P_g,rho_g_oxides,levels=density_levels,alpha=1,cmap=density_cmap)
axi.contour(T_g-273.15,P_g,rho_g_oxides,levels=highlight_densities,alpha=1,colors="black")
plt.xlim(T_limits)
plt.ylim([Pmin,Pmax])
plt.colorbar(mappable=s, location="left", label="density (10$^3$ kg/m$^3$)")
plt.gca().set_title(composition_to_label(composition))
save_current_fig_as("density")

# Plot phases

fig = plt.figure(figsize=(12,14))
axi = fig.add_subplot(1,1,1)
scs = []
for i,ph in enumerate(uniquestrs):
    scs.append(axi.scatter(T_g[phaseis_g==i]-273.15,P_g[phaseis_g==i],alpha=0.35,label=ph,s=50))
fig.legend(handles=scs,loc="center left")
plot_phase_labels(axi)
plt.xlabel("Temperature (°C)", labelpad=2)
plt.ylabel("Pressure (GPa)")
plt.xlim(T_limits)
plt.xticks(T_ticks)
plt.gca().set_xticklabels(T_tick_labels)
plt.ylim([Pmin,Pmax])
plt.gca().set_title(composition_to_label(composition))
save_current_fig_as("phases")


# Plot variance

fig = plt.figure(figsize=(12,14))
axi = fig.add_subplot(1,1,1)
max_variance = np.nanmax(variance_g)
min_variance = np.nanmin(variance_g)
#s = axi.scatter(T_g-273.15,P_g,c=variance_g,s=50,cmap=variance_cmap,alpha=0.5,vmin=min_variance-1,vmax=max_variance+1)
s=axi.imshow(variance_g,cmap=variance_cmap, vmin=min_variance-1, vmax=max_variance+1,**imshow_kwargs)
plot_phase_labels(axi)
plt.xlabel("Temperature (°C)", labelpad=2)
plt.ylabel("Pressure (GPa)")
plt.xlim(T_limits)
plt.xticks(T_ticks)
plt.gca().set_xticklabels(T_tick_labels)
plt.ylim([Pmin,Pmax])
plt.gca().set_title(composition_to_label(composition))
save_current_fig_as("phase-variance")

# Plot solution time

fig = plt.figure(figsize=(12,14))
axi = fig.add_subplot(1,1,1)
s = plt.imshow(stime_g, cmap=stime_cmap,**imshow_kwargs)
#s = plt.scatter(T_g-273.15,P_g,c=stime_g,s=100,alpha=0.75,cmap=stime_cmap)#,norm=mpl.colors.LogNorm(),)
fig.colorbar(s,location="left",ax=axi, label="sol time (s)")
plt.xlabel("Temperature (°C)", labelpad=2)
plt.ylabel("Pressure (GPa)")
plt.xlim(T_limits)
plt.xticks(T_ticks)
plt.gca().set_xticklabels(T_tick_labels)
plt.ylim([Pmin,Pmax])
plt.gca().set_title(composition_to_label(composition))
save_current_fig_as("stime")

# Plot comparison with Perple_X density

interp = ppx_rho_interpolator(composition,"1000kg/m3")
fig = plt.figure(figsize=(10,5))

# Panel (1): contour reactive density
axi = fig.add_subplot(1,3,1)
s = axi.contourf(T_g-273.15, P_g, rho_g, levels=density_levels, alpha=0.75, cmap=density_cmap)
axi.contour(T_g-273.15, P_g, rho_g, levels=density_levels, alpha=1, cmap=density_cmap)
plt.xlim(T_limits)
plt.xticks(T_ticks)
plt.ylabel("Pressure (GPa)")
plt.gca().set_xticklabels(T_tick_labels)
plt.colorbar(mappable=s, location="bottom", ticks=density_ticks, label="Density (10$^3$ kg/m$^3$)")
plt.gca().set_title(composition_to_label(composition))

if (interp is not None):
    rho_eq_g = interp((T_g, P_g))

    # Panel (2): contour Eqm density
    axi = fig.add_subplot(1,3,2)
    s = axi.contourf(T_g-273.15, P_g, rho_eq_g, levels=density_levels, alpha=0.75, cmap=density_cmap)
    axi.contour(T_g-273.15, P_g, rho_eq_g, levels=density_levels, alpha=1, cmap=density_cmap)
    plt.xlim(T_limits)
    plt.xticks(T_ticks)
    plt.xlabel("Temperature (°C)", labelpad=2)
    plt.gca().set_xticklabels(T_tick_labels)
    plt.ylim([Pmin,Pmax])
    plt.yticks([])
    plt.colorbar(mappable=s, location="bottom",ticks=density_ticks, label="Density (10$^3$ kg/m$^3$)")
    plt.gca().set_title("Equilibrium system (NCFMAS)")

    # Panel (3): Diff b/w Eqm and reactive
    axi = fig.add_subplot(1,3,3)
    diff = (rho_g-rho_eq_g)/rho_eq_g*100
    absmax = np.ceil(np.nanmax(np.absolute(diff)))
    levels = np.arange(-absmax, absmax+1, 1)
    s=axi.imshow(diff,cmap=diff_cmap,vmin=-absmax,vmax=absmax, **imshow_kwargs)
    axi.contour(T_g-273.15 ,P_g, diff, levels=levels,**contour_kwargs)
    plt.yticks([])
    plt.colorbar(mappable=s, location="bottom", label="Relative error (%)")
    plt.xticks(T_ticks)
    plt.gca().set_xticklabels(T_tick_labels)
    plt.gca().set_title("Error")

plt.tight_layout()
save_current_fig_as("density_comparison")

# Plot comparison with pyrolite (including oxides)

fig = plt.figure(figsize=(10,5))

axi = fig.add_subplot(1,3,1)
s = axi.contourf(T_g-273.15, P_g, rho_g_oxides, levels=density_levels, alpha=0.75, cmap=density_cmap)
axi.contour(T_g-273.15, P_g, rho_g_oxides, levels=density_levels, alpha=1, cmap=density_cmap)
plt.xlim(T_limits)
plt.xticks(T_ticks)
plt.ylabel("Pressure (GPa)")
plt.gca().set_xticklabels(T_tick_labels)
plt.ylim([Pmin,Pmax])
plt.colorbar(mappable=s, location="bottom", ticks=density_ticks, label="Density (10$^3$ kg/m$^3$)")
plt.gca().set_title(composition_to_label(composition))

axi = fig.add_subplot(1,3,2)
s = axi.contourf(T_g-273.15, P_g, rho_pyrolite_g, levels=density_levels, alpha=0.75, cmap=density_cmap)
axi.contour(T_g-273.15, P_g, rho_pyrolite_g, levels=density_levels, alpha=1, cmap=density_cmap)
plt.xlim(T_limits)
plt.xticks(T_ticks)
plt.xlabel("Temperature (°C)", labelpad=2)
plt.gca().set_xticklabels(T_tick_labels)
plt.ylim([Pmin,Pmax])
plt.yticks([])
plt.colorbar(mappable=s, location="bottom", ticks=density_ticks, label="Density (10$^2$ kg/m$^3$)")
plt.gca().set_title("Pyrolite (Xu et al. 2008)")

axi = fig.add_subplot(1,3,3)
diff = (rho_g_oxides-rho_pyrolite_g)*1000 # kg/m3
maxval = np.ceil(np.nanmax(np.absolute(diff))/100)*100+50
levels = np.arange(-maxval,maxval+50,50)
s=axi.imshow(diff,cmap=diff_cmap,vmin=-maxval,vmax=maxval, **imshow_kwargs)
axi.contour(T_g-273.15, P_g, diff, levels=levels,**contour_kwargs)
plt.colorbar(mappable=s, location="bottom", label="kg/m$^3$")
plt.gca().set_title("Difference")
plt.xticks(T_ticks)
plt.gca().set_xticklabels(T_tick_labels)
plt.yticks([])

plt.tight_layout()
save_current_fig_as("density_contrast")


# Plot combined figure

fig = plt.figure(figsize=(10,5))

# Plot (1,1): density of reactive system (including oxides)
axi = fig.add_subplot(1,3,1)
s = axi.contourf(T_g-273.15, P_g, rho_g_oxides, levels=density_levels, alpha=0.75, cmap=density_cmap)
axi.contour(T_g-273.15, P_g, rho_g_oxides, levels=density_levels, alpha=1, cmap=density_cmap)
plt.xlim(T_limits)
plt.xticks(T_ticks)
plt.gca().set_xticklabels(T_tick_labels)
plt.ylim([Pmin,Pmax])
plt.yticks(P_ticks)
plt.colorbar(mappable=s, location="bottom", ticks=density_ticks, label="Density (10$^3$ kg/m$^3$)",pad=0.12)
plt.ylabel("Pressure (GPa)")
plt.gca().set_title("Reactive system")

# Plot (1,2): difference with equilibrium (doesn't include oxides)
axi = fig.add_subplot(1,3,2)
diff = (rho_g-rho_eq_g)/rho_eq_g*100
absmax = np.ceil(np.nanmax(np.absolute(diff))) # rounded to nearest digit
levels = np.arange(-absmax, absmax+0.2, 0.2)
label_levels = [l for i,l in enumerate(levels) if i%5==0]
s=axi.imshow(diff,cmap=diff_cmap,vmin=-absmax,vmax=absmax, **imshow_kwargs)
contours = axi.contour(T_g-273.15 ,P_g, diff, levels=levels,**contour_kwargs)
axi.clabel(contours, label_levels, inline=True, fmt=rel_err_fmt, fontsize=SMALL_SIZE)
plt.colorbar(mappable=s, location="bottom", label="rel. error (%)",pad=0.12)
plt.xticks(T_ticks)
plt.gca().set_xticklabels(T_tick_labels)
plt.yticks([])
plt.xlabel("Temperature (°C)", labelpad=2)
plt.gca().set_title("Error with equilibrium")

# Plot (1,3): difference with pyrolite (including oxides)
axi = fig.add_subplot(1,3,3)
diff = (rho_g_oxides-rho_pyrolite_g)*1000 # kg/m3
maxval = np.ceil(np.nanmax(np.absolute(diff))/100)*100+40
levels = np.arange(-maxval,maxval+20,20)
s=axi.imshow(diff,cmap=diff_cmap,vmin=-maxval,vmax=maxval, **imshow_kwargs)
contours = axi.contour(T_g-273.15, P_g, diff, levels=levels,**contour_kwargs)
contour0 = axi.contour(T_g-273.15, P_g, diff, levels=highlight_diffs,**contour_highlight_kwargs)
axi.clabel(contours, contours.levels, inline=True, fmt=diff_fmt, fontsize=SMALL_SIZE)
axi.clabel(contour0, contour0.levels, inline=True, fmt=diff_fmt, fontsize=SMALL_SIZE)
plt.colorbar(mappable=s, location="bottom", label="Difference (kg/m$^3$)", pad=0.12)
plt.gca().set_title("Density above pyrolite")
plt.xticks(T_ticks)
plt.gca().set_xticklabels(T_tick_labels)
plt.yticks([])

plt.tight_layout()

save_current_fig_as("density_summary")