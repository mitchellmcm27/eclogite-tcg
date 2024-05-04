import numpy as np
from contextlib import redirect_stdout, redirect_stderr
from thermoengine.equilibrate import System
from thermoengine.samples import Assemblage, SampleMaker
from thermoengine.chemistry import OxideMolComp
import time
import io

from ..base import *

import matplotlib.pyplot as plt

class EquilibratePD:

    T    = None
    p    = None
    comp = None

    system = None
    assemblage = None
    phases = None
    options = None
    
    stdout = ''
    stderr = ''
    excstr = ''
    stime  = None
    
    def __init__(self, phases):
        self.phases = phases
        
    def solve(self,T,p,comp,**kwargs):
        self.set_initial_params(T,p,comp,**kwargs)

        try:
          so = io.StringIO()
          se = io.StringIO()
          with redirect_stdout(so), redirect_stderr(se):
            tic = time.perf_counter()
            if self.system is not None and 'options' not in kwargs:
                self.system.update(T=T, P=p, comp=comp)
            else:
                self.system = System(T=T, P=p, comp=comp,
                                     options=self.options,
                                     phase_library=self.phases)
            self.assemblage = self.system.stable_assemblage
            toc = time.perf_counter()
            self.stime = toc-tic
          self.stdout = so.getvalue()
          self.stderr = se.getvalue()
        except Exception as e:
          self.assemblage  = None
          self.excstr = repr(e)

    def set_initial_params(self,T,p,comp,**kwargs):
        self.T     = T
        self.p     = p
        self.comp  = comp
        self.stdout = kwargs.get('stdout', '')
        self.stderr = kwargs.get('stderr', '')
        self.excstr = kwargs.get('excstr', '')
        self.stime  = kwargs.get('stime', None)
        self.assemblage  = kwargs.get('assemblage', None)
        self.options = kwargs.get('options', {'grid_spacing':1/10})

    def final_phases(self):
        assert(self.assemblage is not None)
        phaseabbrev = sorted(self.assemblage.sample_names.tolist())
        phasenames = [self.phases.available_phase_names[self.phases.available_phase_abbrevs.index(abbrev)] \
                                                          for abbrev in phaseabbrev]
        return phasenames, phaseabbrev
    
    def final_rho(self):
        assert(self.assemblage is not None)

        phaseabbrev = self.assemblage.sample_names
        phaseobjs = [self.phases._available_phases[self.phases._available_phase_abbrevs.index(abbrev)] for abbrev in phaseabbrev]

        endmemcomps = self.assemblage.sample_endmem_comps
        mols = [emcomp if len(emcomp) > 1 else None for emcomp in endmemcomps]
        volumes = np.asarray([ph.volume(self.T, self.p, mol=mols[i]) for i, ph in enumerate(phaseobjs)])
        masses = np.asarray([np.dot(ph.props['molwt'], endmemcomps[i]) for i, ph in enumerate(phaseobjs)])
        rhois = masses/volumes

        mols = self.assemblage.sample_amounts
        totmass = np.dot(mols, masses)
        massfracs = np.asarray([m*mols[i]/totmass for i, m in enumerate(masses)])

        rho = 1./sum(massfracs/rhois)
        return rho
    
    def final_Ci0(self):
        assert(self.assemblage is not None)

        phaseabbrev = self.assemblage.sample_names
        phaseobjs = [self.phases._available_phases[self.phases._available_phase_abbrevs.index(abbrev)] for abbrev in phaseabbrev]

        endmemcomps = self.assemblage.sample_endmem_comps
        masses = np.asarray([np.dot(ph.props['molwt'], endmemcomps[i]) for i, ph in enumerate(phaseobjs)])
        C = [[ph.props['molwt'][k]*endmemcomps[i][k]/masses[i] for k in range(len(endmemcomps[i]))] for i, ph in enumerate(phaseobjs)]

        return [Ci[0] for Ci in C]
    
    def final_Xi0(self):
        assert(self.assemblage is not None)

        X = self.assemblage.sample_endmem_comps

        return [Xi[0] for Xi in X]
    
#########################################################################

class EquilibratePDGrid:
    """A class that contains a grid of reactive phase diagram equilibrate objects.
    
    It has been specially constructed to avoid saving the assemblage object, which cannot be pickled, hence a lot of output is stored in nested lists."""
    ax_args = None

    x_range = None
    y_range = None

    xgrid   = None
    ygrid   = None

    Tgrid    = None
    pgrid    = None
    cgrid = None

    equilgrid = None
    
    kwargs = None
    
    def solve(self,phases,ax_args,x_range,y_range,**kwargs):
        self.ax_args = ax_args

        self.x_range = x_range
        self.y_range = y_range

        self.kwargs = kwargs

        T = kwargs.pop('T', None)
        p = kwargs.pop('p', None)
        comp = kwargs.pop('comp', None)
        
        self.xgrid, self.ygrid = np.meshgrid(self.x_range, self.y_range)
        
        self.Tgrid     = [[None for j in range(len(self.x_range))] for i in range(len(self.y_range))]
        self.pgrid     = [[None for j in range(len(self.x_range))] for i in range(len(self.y_range))]
        self.cgrid     = [[None for j in range(len(self.x_range))] for i in range(len(self.y_range))]

        self.phnamegrid = [[[] for j in range(len(self.x_range))] for i in range(len(self.y_range))]
        self.amountgrid = [[[] for j in range(len(self.x_range))] for i in range(len(self.y_range))]
        self.emcompgrid = [[[] for j in range(len(self.x_range))] for i in range(len(self.y_range))]

        self.outgrid   = [[''   for j in range(len(self.x_range))] for i in range(len(self.y_range))]
        self.errgrid   = [[''   for j in range(len(self.x_range))] for i in range(len(self.y_range))]
        self.excgrid   = [[''   for j in range(len(self.x_range))] for i in range(len(self.y_range))]
        self.stimegrid = [[None for j in range(len(self.x_range))] for i in range(len(self.y_range))]

        Xmin = 0.0
        Xmax = 1.0
        if 'X' in ax_args:
            ax = ax_args.index('X')
            if ax==0:
                X_range = x_range
            else:
                X_range = y_range
            Xmin = X_range[0]
            Xmax = X_range[-1]
        
        for i,y in enumerate(self.y_range):
            for j,x in enumerate(self.x_range):
                equil = EquilibratePD(phases)

                vals = [x,y]
                for ax,arg in enumerate(ax_args):
                    if arg == 'T': 
                        T = vals[ax]
                    if arg == 'p':
                        p = vals[ax]
                    if arg == 'X':
                        X0 = kwargs['X0']
                        X1 = kwargs['X1']
                        f = vals[ax]
                        assert(f >= 0.0)
                        assert(f <= 1.0)
                        comp = OxideMolComp()
                        cdict = {k : X0.get(k, 0.0)*f + X1.get(k, 0.0)*(1.0-f) for k in comp.all_components}
                        comp = OxideMolComp(**cdict)
                
                equil.solve(T,GPa2Bar(p),comp,**kwargs)
                
                self.Tgrid[i][j]    = T
                self.pgrid[i][j]    = p
                self.cgrid[i][j]    = comp.all_data

                if equil.assemblage is not None:
                    self.phnamegrid[i][j] = equil.assemblage.sample_names
                    self.amountgrid[i][j] = equil.assemblage.sample_amounts
                    self.emcompgrid[i][j] = equil.assemblage.sample_endmem_comps

                self.outgrid[i][j] = equil.stdout
                self.errgrid[i][j] = equil.stderr
                self.excgrid[i][j] = equil.excstr
                self.stimegrid[i][j] = equil.stime

class PDEquilibrateGridDiagnostics:
    phasenames = None
    phasestrs = None
    uniquestrs = None
    phaseis = None
    
    def __init__(self, phases, grid):
        self.phases = phases
        self.grid = grid
        
    def reconstruct_equil(self,i,j):
        T      = self.grid.Tgrid[i][j]
        p      = self.grid.pgrid[i][j]
        cdict  = self.grid.cgrid[i][j]

        phname = self.grid.phnamegrid[i][j]
        amount = self.grid.amountgrid[i][j]
        emcomp = self.grid.emcompgrid[i][j]

        phaseobjs = [self.phases._available_phases[self.phases._available_phase_abbrevs.index(abbrev)] for abbrev in phname]
        
        assemblage = Assemblage([SampleMaker.get_fixed_comp_sample(ph, X=emcomp[pi], T=T, P=p, amount=amount[pi]) for pi, ph in enumerate(phaseobjs)])

        stdout = self.grid.outgrid[i][j]
        stderr = self.grid.errgrid[i][j]
        excstr = self.grid.excgrid[i][j]
        stime  = self.grid.stimegrid[i][j]

        equil = EquilibratePD(self.phases)
        equil.set_initial_params(T,p,cdict,assemblage=assemblage,stdout=stdout,stderr=stderr,excstr=excstr,stime=stime)

        return equil

    def phase_diagnostics(self):
        self.phasenames = [[[] for j in range(len(self.grid.x_range))] for i in range(len(self.grid.y_range))]
        self.phasestrs = [['' for j in range(len(self.grid.x_range))] for i in range(len(self.grid.y_range))]
        for i,P in enumerate(self.grid.y_range):
            for j,x in enumerate(self.grid.x_range):
                equil = self.reconstruct_equil(i,j)
                if equil.assemblage is not None:
                    equilphasenames, phaseabbrev = equil.final_phases()
                    self.phasenames[i][j] = equilphasenames
                    self.phasestrs[i][j] = '+'.join(phaseabbrev)
        
        # a set of unique strings representing the phase combinations present across the phase diagram
        self.uniquestrs = sorted(list(set([phstr for phstrl in self.phasestrs for phstr in phstrl if phstr != ''])))
        
        def index(ls,v):
            try:
                i = ls.index(v)
            except ValueError:
                i = -1
            return i


        # an index into the uniquestr list to give a grid showing which phases are present where
        self.phaseis = np.asarray([[index(self.uniquestrs,phstr) for phstr in phasestrr]  for phasestrr in self.phasestrs])

    def update(func):
        def update_decorator(self, *args, **kwargs):
            if self.phaseis is None: self.phase_diagnostics()
            return func(self, *args, **kwargs)
        return update_decorator

    @update
    def scatter(self,ax,x,y,c,cmap=plt.get_cmap('Paired'),norm=None):
        s = ax.scatter(x,y,c=c,s=100,alpha=0.75,cmap=cmap,norm=norm)
        return s
    
    @update
    def plot_phase_labels(self, ax):
        fxgrid = self.grid.xgrid.flatten()
        fygrid = self.grid.ygrid.flatten()
        fphaseis = self.phaseis.flatten()

        for i,ph in enumerate(self.uniquestrs):
            if i in fphaseis:
                indices = fphaseis==i
                ax.text(np.mean(fxgrid[indices]), np.mean(fygrid[indices]),ph,\
                        backgroundcolor='white',fontsize=12)

        
    @update
    def plot_rho(self):
        
        rhogrid = np.empty(self.grid.ygrid.shape)
        for i,P in enumerate(self.grid.y_range):
            for j,x in enumerate(self.grid.x_range):
                equil = self.reconstruct_equil(i,j)
                if equil.assemblage is not None:
                    rhogrid[i][j] = equil.final_rho()
        
        fig = plt.figure(figsize=(12,14))
        axi = fig.add_subplot(1,1,1)
        ax = self.setup_axes(axi)

        cmap = plt.get_cmap('bwr')
        s = self.scatter(ax,self.grid.xgrid[self.phaseis>=0],self.grid.ygrid[self.phaseis>=0],rhogrid[self.phaseis>=0],cmap)
        fig.colorbar(s,location='left',ax=axi)
        
        self.plot_phase_labels(ax)

    @update
    def plot_stime(self):
        fig = plt.figure(figsize=(12,14))
        axi = fig.add_subplot(1,1,1)
        ax = self.setup_axes(axi)

        cmap = plt.get_cmap('bwr')
        s = self.scatter(ax,self.grid.xgrid[self.phaseis>=0],self.grid.ygrid[self.phaseis>=0],np.asarray(self.grid.stimegrid)[self.phaseis>=0],cmap)#,norm=mpl.colors.LogNorm())
        fig.colorbar(s,location='left',ax=axi)

        self.plot_phase_labels(ax)

    def setup_axes(self,axi):
        return axi
    
    @update
    def plot_phases(self):
        fig = plt.figure(figsize=(24,14))
        axes = fig.subplot_mosaic([['A',[['1','2'],['3','4'],['5','6']]]])
        axi = axes['A']
        axi.set_navigate(False)
        ax = self.setup_axes(axi)
        if axi != ax: fig.add_axes(ax)
        ax.set_navigate(True)
        
        raxes = [axes[repr(i)] for i in range(1,7)]
        for axis in raxes: axis.set_visible(False)
    
        jgrid, igrid = np.meshgrid(range(len(self.grid.x_range)), range(len(self.grid.y_range)))
        scs = []
        sis = []
        sjs = []
        for i,ph in enumerate(self.uniquestrs):
            sis.append(igrid[self.phaseis==i])
            sjs.append(jgrid[self.phaseis==i])
            scs.append(ax.scatter(self.grid.xgrid[self.phaseis==i],self.grid.ygrid[self.phaseis==i],alpha=0.5,label=ph,s=50))
        fig.legend(handles=scs,loc='center left')
        
        #self.plot_phase_labels(ax)

        annot = ax.annotate("", xy=(0,0), xytext=(5,5), textcoords='offset points', backgroundcolor='white')
        annot.set_visible(False)
        
        def hover(event):
            if event.inaxes==ax:
                labels = [self.uniquestrs[si]+(','.join([' (i,j)=({},{})'.format(sis[si][i],sjs[si][i]) for i in contains[1]['ind']])) for si, sc in enumerate(scs) if (contains:=sc.contains(event))[0]]
                if len(labels)>0:
                    annot.xy = (event.xdata, event.ydata)
                    annot.set_text('; '.join(labels))
                    annot.set_visible(True)
                else:
                    annot.set_visible(False)               
            else:
                annot.set_visible(False)
        
        fig.canvas.mpl_connect("motion_notify_event", hover)

#    @update
#    def plot_Xi1(self):
#        I = len(self.rxn.phases())
#        Xi1grid = np.empty(self.grid.ygrid.shape+(I,))
#        for i,P in enumerate(self.grid.y_range):
#            for j,x in enumerate(self.grid.x_range):
#                ode = self.reconstruct_ode(i,j)
#                if ode.sol is not None:
#                    Xi1grid[i][j] = np.ones(ode.I)-np.asarray(ode.final_Xi0())
#
#        fig = plt.figure(figsize=(21,np.ceil(I/3)*10))
#
#        for i in range(I):
#            name = self.rxn.phases()[i].name()
#            phasein = np.asarray([[name in self.phasenames[i][j] for j in range(len(self.grid.x_range))] for i in range(len(self.grid.y_range))])
#
#            axi = fig.add_subplot(int(np.ceil(I/3)),3,i+1)
#            ax = self.setup_axes(axi)
#            ax.set_title(name)
#
#            cmap = plt.get_cmap('bwr')
#            s = self.scatter(ax,self.grid.xgrid[phasein],self.grid.ygrid[phasein],Xi1grid[phasein][:,i],cmap=cmap)
#            fig.colorbar(s,location='left',ax=axi)
#            
#    @update
#    def plot_Ci1(self):
#        I = len(self.rxn.phases())
#        Ci1grid = np.empty(self.grid.ygrid.shape+(I,))
#        for i,P in enumerate(self.grid.y_range):
#            for j,x in enumerate(self.grid.x_range):
#                ode = self.reconstruct_ode(i,j)
#                if ode.sol is not None:
#                    Ci1grid[i][j] = np.ones(ode.I)-np.asarray(ode.final_Ci0())
#
#        fig = plt.figure(figsize=(21,np.ceil(I/3)*10))
#
#        for i in range(I):
#            name = self.rxn.phases()[i].name()
#            phasein = np.asarray([[name in self.phasenames[i][j] for j in range(len(self.grid.x_range))] for i in range(len(self.grid.y_range))])
#
#            axi = fig.add_subplot(int(np.ceil(I/3)),3,i+1)
#            ax = self.setup_axes(axi)
#            ax.set_title(name)
#
#            cmap = plt.get_cmap('bwr')
#            s = self.scatter(ax,self.grid.xgrid[phasein],self.grid.ygrid[phasein],Ci1grid[phasein][:,i],cmap=cmap)
#            fig.colorbar(s,location='left',ax=axi)
    
