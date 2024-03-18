import io
import numpy as np
from scipy.integrate import solve_ivp
from contextlib import redirect_stdout, redirect_stderr

from ..base import *

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib._color_data as mcd
colors = list(mcd.TABLEAU_COLORS.values())+['#000000']
acolors = ['c','m','y','k']
from cycler import cycler
dcycler = (cycler(linestyle=['-', '--', ':', '-.'])*cycler(color=colors))
acycler = (cycler(linestyle=['-', '--', ':', '-.'])*cycler(color=['c', 'm', 'y', 'k']))

class BasePDReactiveODE:
    rxn  = None

    I = 0
    Kis = []
    K = 0
    J = 0
    
    T    = None
    p    = None
    Da   = None
    rho0 = None
    eps  = None
    F0   = None
    Cik0 = None
    
    sol  = None
    stdout = ''
    stderr = ''
    excstr = ''
    stime  = None
    
    def __init__(self,rxn):
 
        self.rxn = rxn
        self.I = len(self.rxn.phases())
        self.Kis = [len(self.rxn.phases()[i].endmembers()) for i in range(self.I)]
        self.K = sum(self.Kis)
        self.J = len(self.rxn.nu())
        
    def set_initial_params(self,T,p,mi0,Cik0,**kwargs):
        self.eps   = kwargs.get('eps', 1.e-2)
        self.Da    = kwargs.get('Da', 1.0)

        self.T     = T
        self.p     = p
        self.mi0   = mi0
        self.Cik0  = Cik0
        C  = self.reshapeC(self.Cik0)
        Cs = self.regularizeC(C)
        self.rho0  = 1./self.v(self.T, self.p, Cs, self.mi0)

        self.sol    = kwargs.get('sol', None)
        self.stdout = kwargs.get('stdout', '')
        self.stderr = kwargs.get('stderr', '')
        self.excstr = kwargs.get('excstr', '')
        self.stime  = kwargs.get('stime', None)

    def final_phases(self,tol):
        assert(self.sol is not None)
        phasenames = np.array([self.rxn.phases()[i].name() for i in range(self.I)], dtype=str)
        phaseabbrev = np.array(['Pv' if (abbrev := self.rxn.phases()[i].abbrev()) == "MgFePv" else abbrev for i in range(self.I)], dtype=str)
        mi = self.sol.y[:self.I,-1]
        phaseindices = mi > tol
        return phasenames[phaseindices], phaseabbrev[phaseindices]
    
    def final_rho(self):
        assert(self.sol is not None)
        mi  = self.sol.y[:self.I,-1]
        Cik = self.sol.y[self.I:self.I+self.K,-1]

        C  = self.reshapeC(Cik)
        Cs = self.regularizeC(C)

        v = self.v(self.T, self.p, Cs, mi)
        rho   = 1/v
        return rho
    
    def final_Ci0(self):
        assert(self.sol is not None)
        Cik = self.sol.y[self.I:self.I+self.K,-1]
        
        C = self.reshapeC(Cik)
        return [Ci[0] for Ci in C]
    
    def final_Xi0(self):
        assert(self.sol is not None)
        Cik = self.sol.y[self.I:self.I+self.K,-1]
        
        C = self.reshapeC(Cik)
        X = self.rxn.C_to_X(C)
        return [Xi[0] for Xi in X]
    
    def reshapeC(self,Cik):
        C = self.rxn.zero_C()
        Kis = 0
        for i,Ki in enumerate(self.Kis):
            C[i] = Cik[Kis:Kis+Ki]
            Kis = Kis+Ki
        return C
    
    def regularizem(self, mi):
        mi = np.asarray(mi)
        mi = mi + self.eps
        return mi
    
    def regularizeC(self,C):
        C = [np.maximum(np.asarray(C[i]), self.eps*np.ones(len(C[i]))) for i in range(len(C))]
        C = [np.asarray(C[i])/sum(C[i]) for i in range(len(C))]
        return C
    
    def rhoi(self, T, p, Cs):
#   void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x, const ufc::cell &cell) const
#   { 
#     update(x, cell);
#     //update_mi(x, cell);
#     zero(values);
    
#     //Evaluate
#     rxn.rho(T,P,cik,rhoi);
#     for (int i=0; i < I; i++)
#     { 
#       values[i] = rhoi[i];
#     }
#   }
        rhoi = np.asarray(self.rxn.rho(T, p, Cs))
        return rhoi
        
    def drhoidus(self, T, p, Cs):
#   void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x, const ufc::cell &cell) const
#   {
#     update(x, cell);
#     //update_mi(x, cell);
#     zero(values);

#     //Evaluate
#     int sKi = 0;
#     for (int i=0; i < I; i++)
#     {
#       std::vector<double> drhoidcik = rxn.phases()[i]->drho_dc(T,P,cik[i]);
#       for (int di=0; di < I; di++)
#       {
#         values[i*(I+K)+di] = 0;
#       }
#       for (int dk=0; dk < Kis[i]; dk++)
#       {
#         values[i*(I+K) + I + sKi + dk] = drhoidcik[dk];
#       }
#       sKi += Kis[i];
#     }
#   }
        drhoidus = np.zeros((self.I,self.I+self.K))
        sKi = 0
        for i in range(self.I):
            drhoidcik = self.rxn.phases()[i].drho_dc(T, p, Cs[i])
            for dk in range(self.Kis[i]):
                drhoidus[i, self.I + sKi + dk] = drhoidcik[dk]
            sKi += self.Kis[i]
#         print('drhoidus:')
#         print('')
#         for i,v in enumerate(drhoidus.flatten()):
#             print('{}: {}'.format(i,v))
        return drhoidus
    
    def v(self, T, p, Cs, mi):
        rhoi = self.rhoi(T,p,Cs)
        v    = sum(mi/rhoi)
        return v
    
    def dvdus(self, T, p, Cs, mi):
        rhoi = self.rhoi(T,p,Cs)
        drhoidus = self.drhoidus(T,p,Cs)
        dvdus = np.zeros(self.I+self.K)
        sKi = self.I
        for di in range(self.I):
            dvdus[di] = 1./rhoi[di]
            for dk in range(self.Kis[di]):
                dvdus[sKi+dk] = -sum(mi/(rhoi**2)*drhoidus[:,sKi+dk])
            sKi += self.Kis[di]
        return dvdus
    
    def Gammai(self, T, p, Cs, mi):
#   void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x, const ufc::cell &cell) const
#   {
#     update(x, cell);
#     update_mi(x, cell);
#     zero(values);

#     //Evaluate
#     rxn.Gamma_i(T,P,cik,mi,Gammai);
#     //for (int i=0; i < I-1; i++)
#     for (int i=0; i < I; i++)
#     {
#       values[i] = Gammai[i];
#     }
#   }
        gammai  = np.asarray(self.rxn.Gamma_i(T,p,Cs,mi))
        return gammai
    
    def dGammaidus(self, T, p, Cs, mi):
#   void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x, const ufc::cell &cell) const
#   { 
#     update(x, cell);
#     update_mi(x, cell);
#     zero(values);
    
#     //Evaluate
#     rxn.dGamma_i_dPhi(T,P,cik,mi,dGammaidmi);
#     rxn.dGamma_i_dC(T,P,cik,mi,dGammaidcik);
#     //for (int i=0; i < I-1; i++)
#     int sKi = 0;
#     for (int i=0; i < I; i++)
#     { 
#       for (int di=0; di < I; di++)
#       { 
#         values[sKi+di] = dGammaidmi[i][di];
#       }
#       sKi += I;
#       for (int di=0; di < I; di++)
#       { 
#         for (int dk=0; dk < Kis[di]; dk++)
#         { 
#           values[sKi+dk] = dGammaidcik[i][di][dk];
#         }
#         sKi += Kis[di];
#       }
#     }
#   }        
        dGammaidmi = np.asarray(self.rxn.dGamma_i_dPhi(T,p,Cs,mi))
        dGammaidCik = self.rxn.dGamma_i_dC(T,p,Cs,mi)
        dGammaidus = np.zeros((self.I, self.I+self.K))
        for i in range(self.I):
            for di in range(self.I):
                dGammaidus[i][di] = dGammaidmi[i][di]
            sKi = self.I
            for di in range(self.I):
                for dk in range(self.Kis[di]):
                    dGammaidus[i][sKi+dk] = dGammaidCik[i][di][dk]
                sKi += self.Kis[di]
#         print('dGammaidus:')
#         print('')
#         for i,v in enumerate(dGammaidus.flatten()):
#             print('{}: {}'.format(i,v))
        return dGammaidus
    
    def Gammaik(self, T, p, Cs, mi):
#   void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x, const ufc::cell &cell) const
#   { 
#     update(x, cell);
#     update_mi(x, cell);
#     zero(values);
    
#     //Evaluate
#     rxn.Gamma_ik(T,P,cik,mi,Gammaik);
#     //for (int i=0; i < I-1; i++)
#     int sKi = 0;
#     for (int i=0; i < I; i++)
#     { 
#       for (int k=0; k < Kis[i]; k++)
#       { 
#         values[sKi+k] = Gammaik[i][k];
#       }
#       sKi += Kis[i];
#     }
#   }
        gamma_ik = self.rxn.Gamma_ik(T,p,Cs,mi)
        Gammaik = np.zeros(self.K)
        sKi = 0
        for i in range(self.I):
            for k in range(self.Kis[i]):
                Gammaik[sKi+k] = gamma_ik[i][k]
            sKi += self.Kis[i]
        return Gammaik
    
    def dGammaikdus(self, T, p, Cs, mi):
#   void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x, const ufc::cell &cell) const
#   {
#     update(x, cell);
#     update_mi(x, cell);
#     zero(values);

#     //Evaluate
#     //for (int i=0; i < I-1; i++)
#     int sKi = 0;
#     for (int i=0; i < I; i++)
#     {
#       rxn.dGamma_ik_dPhi(T,P,cik,mi,i,dGammaikdmii);
#       for (int k=0; k < Kis[i]; k++)
#       {
#         for (int di=0; di < I; di++)
#         {
#           values[sKi+di] = dGammaikdmii[k][di];
#         }
#         sKi += I;
#         for (int di=0; di < I; di++)
#         {
#           for (int dk=0; dk < Kis[di]; dk++)
#           {
#             dGammaikdcikik = rxn.dGamma_ik_dC(T,P,cik,mi,i,k,di,dk);
#             values[sKi+dk] = dGammaikdcikik;
#           }
#           sKi += Kis[di];
#         }
#       }
#     }
#   }
        dGammaikdus = np.zeros((self.K,self.I+self.K))
        sKi = 0
        for i in range(self.I):
            dGammaikdmii = self.rxn.dGamma_ik_dPhi(T,p,Cs,mi,i)
            for k in range(self.Kis[i]):
                for di in range(self.I):
                    dGammaikdus[sKi+k,di] = dGammaikdmii[k][di]
                sKi2 = self.I
                for di in range(self.I):
                    for dk in range(self.Kis[di]):
                        dGammaikdCikik = self.rxn.dGamma_ik_dC(T,p,Cs,mi,i,k,di,dk)
                        dGammaikdus[sKi+k,sKi2+dk] = dGammaikdCikik
                    sKi2 += self.Kis[di]
            sKi += self.Kis[i]
#         print('dGammaikdus:')
#         print('')
#         for i,v in enumerate(dGammaikdus.flatten()):
#             print('{}: {}'.format(i,v))
        return dGammaikdus

    def rhs(self, t, u):
        '''
        return rhs of the dimensionaless 1-D reactive equations

        Parameters
        ----------

        t: float
            parameterized distance
        u: array
            solution array [ mi, Cik ]

        '''
#         print('  In rhs')
        # Extract variables
        mi  = u[:self.I]
        Cik = u[self.I:self.I+self.K]

        C  = self.reshapeC(Cik)
        Cs = self.regularizeC(C)
        mis = self.regularizem(mi)
        T = self.T
        p = self.p
        
        v = self.v(T,p,Cs,mi)
        Gammai = self.Gammai(T,p,Cs,mi)
        Gammaik = self.Gammaik(T,p,Cs,mi)
        
        du = np.zeros(self.I+self.K)
        sKi = 0
        for i in range(self.I):
            du[i] = self.Da*self.rho0*Gammai[i]*v
            for k in range(self.Kis[i]):
                GikcGi = Gammaik[sKi+k] - C[i][k]*Gammai[i]
                du[self.I+sKi+k] = self.Da*self.rho0*GikcGi*v/mis[i]
            sKi += self.Kis[i]
        
        return du
    
    def jac(self, t, u):
        '''
        return rhs of the dimensionaless 1-D reactive equations

        Parameters
        ----------

        t: float
            parameterized distance
        u: array
            solution array [ mi, Cik ]

        '''
#         print('  In jac')
        # Extract variables
        mi  = u[:self.I]
        Cik = u[self.I:self.I+self.K]

        C  = self.reshapeC(Cik)
        Cs = self.regularizeC(C)
        mis = self.regularizem(mi)
        T = self.T
        p = self.p
        
        v = self.v(T,p,Cs,mi)
        dvdus = self.dvdus(T,p,Cs,mi)
        Gammai = self.Gammai(T,p,Cs,mi)
        Gammaik = self.Gammaik(T,p,Cs,mi)
        dGammaidus = self.dGammaidus(T,p,Cs,mi)
        dGammaikdus = self.dGammaikdus(T,p,Cs,mi)
        
        dudu = np.zeros((self.I+self.K,self.I+self.K))
        sKi = 0
        for i in range(self.I):
            dudu[i,:] = self.Da*self.rho0*(dGammaidus[i]*v \
                                         + Gammai[i]*dvdus)
            for k in range(self.Kis[i]):
                GikcGi = Gammaik[sKi+k] - C[i][k]*Gammai[i]
                dGikcGi = dGammaikdus[sKi+k] - C[i][k]*dGammaidus[i]
                dGikcGi[self.I+sKi+k] -= Gammai[i]
                dudu[self.I+sKi+k,:] = self.Da*self.rho0*(dGikcGi*v/mis[i] \
                                                        + GikcGi*dvdus/mis[i])
                dudu[self.I+sKi+k,i] -= self.Da*self.rho0*(GikcGi*v/(mis[i]**2))
            sKi += self.Kis[i]
        
        return dudu
    
    def plot(self):
        fig = plt.figure(figsize=(12,24))
        axes = [fig.add_subplot(3,2,i+1) for i in  range(6)]
        
        self.plot_axes(axes)
        self.apply_labels(axes)

        fig.suptitle('$p$ = {:.1f} GPa, $T$ = {:.1f} K, $\\rho_0$ = {:.1f} kg/m$^3$, eps = {:.2e}'.format(Bar2GPa(self.p), self.T, rho2kgpm3(self.rho0),self.eps),y=0.9)
        
        plt.show()
    
    def apply_labels(self,axes):
        
        plabels = [self.rxn.phases()[i].abbrev() for i in range(self.I)]
        elabels = [self.rxn.phases()[i].endmembers()[k].formula()+'_('+self.rxn.phases()[i].abbrev()+')' \
                                           for i in range(self.I) for k in range(self.Kis[i])]

        axes[0].set_ylabel('$\\phi_i$')
        labels = [line.get_label()  for line in  axes[0].get_lines()[len(plabels):]]
        axes[0].legend(plabels+labels,loc='best')
#         axes[0].set_ylim([-0.05,1.05])
        
        axes[1].set_ylabel('$\Gamma_i$')
        
        axes[2].set_ylabel('$m_i$')
        
        axes[3].set_ylabel('$\\rho$ kg/m$^3$')
        axes[3].legend(loc='best')
        
        axes[4].set_ylabel('$C_i^k$')
        axes[4].set_ylim([-0.05,1.05])
        axes[4].legend(elabels, loc='best')
        
        axes[5].set_ylabel('$\Gamma_i^k$')
        
        for axis in axes:
            axis.grid()
            axis.set_xlabel('$t^*$')
                
    def plot_axes(self, axes, ls='-', label=None):
        assert(self.sol is not None)
        assert(len(axes)==6)
        
        sol = self.sol
        
        t = self.sol.t
        
        T = self.T
        p = self.p
        
        mi  = self.sol.y[:self.I].T
        Cik = self.sol.y[self.I:self.I+self.K].T

        indices_eval = list(range(len(t)))
        gamma_i  = np.zeros(mi.shape)
        gamma_ik = np.zeros(Cik.shape)
        rhoi     = gamma_i.copy()
        rho      = np.zeros(t.shape)
        C_ik     = np.zeros(Cik.shape)

        for i in range(len(t)):
            C = self.reshapeC(Cik[i])
            Cs = self.regularizeC(C)
            gamma_i[i,:]  = self.rxn.Gamma_i(T,p,Cs,mi[i])
            gamma_ik[i,:] = [gik for gi in self.rxn.Gamma_ik(T,p,Cs,mi[i]) for gik in gi]
            rhoi[i]       = self.rxn.rho(T, p, Cs)
            v             = mi[i]/rhoi[i]
            rho[i]        = 1./v.sum()
            C_ik[i,:]     = [cik for ci in Cs for cik in ci]

        #for i in range(self.I): axes[0].plot(t,rho*mi[:,i]/rhoi[:,i],color=colors[i],ls=ls)
        for i in range(self.I): axes[0].plot(t,rho*mi[:,i]/rhoi[:,i],ls=ls)
        mlabel='$\\sum_i\\phi_i$'
        if label is not None: mlabel=mlabel+' '+label
        sumphiline = axes[0].plot(t,rho*(mi/rhoi).sum(axis=-1),color='k',ls='-',label=mlabel)
        
        #for i in range(self.I): axes[1].plot(t,gamma_i[:,i],color=colors[i],ls=ls)
        for i in range(self.I): axes[1].plot(t,gamma_i[:,i],ls=ls)

#         for i in range(self.I): 
#             for k in range(self.Kis[i]):
#                 axes[4].plot(t,mi[:,i]*Cik[:,sum(self.Kis[:i])+k],ls=ls,color=colors[sum(self.Kis[:i])+k])

        #for i in range(self.I): axes[2].plot(t,mi[:,i],color=colors[i],ls=ls)
        for i in range(self.I): axes[2].plot(t,mi[:,i],ls=ls)

        #for i in range(self.I): axes[3].plot(t,rho2kgpm3(rhoi[:,i]),color=colors[i],ls=ls,label='_')
        for i in range(self.I): axes[3].plot(t,rho2kgpm3(rhoi[:,i]),ls=ls,label='_')
        mlabel='$\\rho$'
        if label is not None: mlabel=mlabel+' '+label
        rholine = axes[3].plot(t,rho2kgpm3(rho),'k',ls=ls,label=mlabel)

        #for ik in range(self.K): axes[4].plot(t,Cik[:,ik],ls=ls,color=colors[len(colors)%(ik+1)])
        for ik in range(self.K): axes[4].plot(t,Cik[:,ik],ls=ls)

        #for ik in range(self.K): axes[5].plot(t,gamma_ik[:,ik],ls=ls,color=colors[len(colors)%(ik+1)])
        for ik in range(self.K): axes[5].plot(t,gamma_ik[:,ik],ls=ls)


#########################################################################

class PDReactiveGrid:
    """A class that contains a grid of reactive phase diagram odes.
    
    It has been specially constructed to avoid saving the reaction object, which cannot be pickled, hence a lot of output is stored in nested lists."""
    odecls = None

    i0  = None

    ax_args = None

    x_range = None
    y_range = None

    end = None
    
    xgrid   = None
    ygrid   = None

    Tgrid    = None
    pgrid    = None
    mi0grid  = None
    Cik0grid = None

    solgrid = None
    outgrid = None
    errgrid = None
    excgrid = None
    stimegrid = None
    
    kwargs = None
    
    def solve(self,rxn,odecls,i0,ax_args,x_range,y_range,end,**kwargs):
        self.odecls = odecls

        self.i0 = i0
        
        self.ax_args

        self.x_range = x_range
        self.y_range = y_range

        self.end = end
        
        self.kwargs = kwargs

        T = kwargs.pop('T', None)
        p = kwargs.pop('p', None)
        Cik0 = kwargs.pop('Cik0', None)
        
        self.xgrid, self.ygrid = np.meshgrid(self.x_range, self.y_range)
        
        self.Tgrid     = [[None for j in range(len(self.x_range))] for i in range(len(self.y_range))]
        self.pgrid     = [[None for j in range(len(self.x_range))] for i in range(len(self.y_range))]
        self.mi0grid   = [[None for j in range(len(self.x_range))] for i in range(len(self.y_range))]
        self.Cik0grid  = [[None for j in range(len(self.x_range))] for i in range(len(self.y_range))]
        self.solgrid   = [[dict for j in range(len(self.x_range))] for i in range(len(self.y_range))]
        self.outgrid   = [[''   for j in range(len(self.x_range))] for i in range(len(self.y_range))]
        self.errgrid   = [[''   for j in range(len(self.x_range))] for i in range(len(self.y_range))]
        self.excgrid   = [[''   for j in range(len(self.x_range))] for i in range(len(self.y_range))]
        self.stimegrid = [[None for j in range(len(self.x_range))] for i in range(len(self.y_range))]
        
        for i,y in enumerate(self.y_range):
            for j,x in enumerate(self.x_range):
                ode = self.odecls(rxn)

                vals = [x,y]
                for ax,arg in enumerate(ax_args):
                    if arg == 'T': 
                        T = vals[ax]
                    if arg == 'p':
                        p = vals[ax]
                    if arg in ['Xi0k0', 'Ci0k0']:
                        if arg=='Xi0k0':
                            phase0 = rxn.phases()[self.i0]
                            Xi0k0 = vals[ax]
                            Xi0 = np.zeros(ode.Kis[self.i0])
                            Xi0[0] = Xi0k0
                            Xi0[1:] = (1.-Xi0k0)/(ode.Kis[self.i0]-1)
                            Ci0k0 = phase0.x_to_c(Xi0)[0]
                        else:
                            Ci0k0 = vals[ax]
                        Cik0 = np.zeros(ode.K)
                        for ip in range(ode.I):
                            if ode.Kis[ip] == 1:
                                 Cik0[sum(ode.Kis[:ip])] = 1.
                            else:
                                 Cik0[sum(ode.Kis[:ip])] = Ci0k0
                                 Cik0[sum(ode.Kis[:ip])+1:sum(ode.Kis[:ip+1])][:] = (1.-Ci0k0)/(ode.Kis[ip]-1.)
                
                # initial conditions
                mi0 = np.zeros(ode.I)
                mi0[self.i0] = 1.0
                
                ode.solve(T,GPa2Bar(p),mi0,Cik0,end,**kwargs)
                
                self.Tgrid[i][j]    = T
                self.pgrid[i][j]    = p
                self.mi0grid[i][j]  = mi0
                self.Cik0grid[i][j] = Cik0

                self.solgrid[i][j] = ode.sol
                self.outgrid[i][j] = ode.stdout
                self.errgrid[i][j] = ode.stderr
                self.excgrid[i][j] = ode.excstr
                self.stimegrid[i][j] = ode.stime

class PDReactiveGridDiagnostics:
    """A class to run diagnostics on a pre-computed phase diagram.

    Separated from the phase diagram calculation to prevent unecessary and slow reevaluations."""
    phasenames = None
    phasestrs = None
    uniquestrs = None
    phaseis = None
    
    def __init__(self, rxn, grid):
        self.rxn  = rxn
        self.grid = grid
        
    def reconstruct_ode(self,i,j):
        T      = self.grid.Tgrid[i][j]
        p      = self.grid.pgrid[i][j]
        mi0    = self.grid.mi0grid[i][j]
        Cik0   = self.grid.Cik0grid[i][j]

        sol    = self.grid.solgrid[i][j]
        stdout = self.grid.outgrid[i][j]
        stderr = self.grid.errgrid[i][j]
        excstr = self.grid.excgrid[i][j]
        stime  = self.grid.stimegrid[i][j]

        ode = self.grid.odecls(self.rxn)
        ode.set_initial_params(T,p,mi0,Cik0,sol=sol,stdout=stdout,stderr=stderr,excstr=excstr,stime=stime)

        return ode

    def phase_diagnostics(self):
        self.phasenames = [[[] for j in range(len(self.grid.x_range))] for i in range(len(self.grid.y_range))]
        self.phasestrs = [['' for j in range(len(self.grid.x_range))] for i in range(len(self.grid.y_range))]
        for i,P in enumerate(self.grid.y_range):
            for j,x in enumerate(self.grid.x_range):
                ode = self.reconstruct_ode(i,j)
                if ode.sol is not None:
                    odephasenames, phaseabbrev = ode.final_phases(1.e-2)
                    self.phasenames[i][j] = odephasenames.tolist()
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
                ode = self.reconstruct_ode(i,j)
                if ode.sol is not None:
                    rhogrid[i][j] = ode.final_rho()
        
        fig = plt.figure(figsize=(12,14))
        axi = fig.add_subplot(1,1,1)
        ax = self.setup_axes(axi)

        cmap = plt.get_cmap('bwr')
        s = self.scatter(ax,self.grid.xgrid[self.phaseis>=0],self.grid.ygrid[self.phaseis>=0],rhogrid[self.phaseis>=0],cmap)
        fig.colorbar(s,location='left',ax=axi)
        
        self.plot_phase_labels(ax)

        
    
    @update
    def plot_mindt(self):
        mindtgrid = np.empty(self.grid.ygrid.shape)
        for i,P in enumerate(self.grid.y_range):
            for j,x in enumerate(self.grid.x_range):
                ode = self.reconstruct_ode(i,j)
                if ode.sol is not None:
                    mindtgrid[i][j] = (ode.sol.t[1:]-ode.sol.t[:-1]).min()
        
        fig = plt.figure(figsize=(12,14))
        axi = fig.add_subplot(1,1,1)
        ax = self.setup_axes(axi)

        cmap = plt.get_cmap('bwr')
        s = self.scatter(ax,self.grid.xgrid[self.phaseis>=0],self.grid.ygrid[self.phaseis>=0],mindtgrid[self.phaseis>=0],cmap,norm=mpl.colors.LogNorm())
        fig.colorbar(s,location='left',ax=axi)

        self.plot_phase_labels(ax)
    
    @update
    def plot_ndt(self):
        ndtgrid = np.empty(self.grid.ygrid.shape)
        for i,P in enumerate(self.grid.y_range):
            for j,x in enumerate(self.grid.x_range):
                ode = self.reconstruct_ode(i,j)
                if ode.sol is not None:
                    ndtgrid[i][j] = len(ode.sol.t)
        
        fig = plt.figure(figsize=(12,14))
        axi = fig.add_subplot(1,1,1)
        ax = self.setup_axes(axi)

        cmap = plt.get_cmap('bwr')
        s = self.scatter(ax,self.grid.xgrid[self.phaseis>=0],self.grid.ygrid[self.phaseis>=0],ndtgrid[self.phaseis>=0],cmap)#,norm=mpl.colors.LogNorm())
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

    @update
    def plot_Xi1(self):
        I = len(self.rxn.phases())
        Xi1grid = np.empty(self.grid.ygrid.shape+(I,))
        for i,P in enumerate(self.grid.y_range):
            for j,x in enumerate(self.grid.x_range):
                ode = self.reconstruct_ode(i,j)
                if ode.sol is not None:
                    Xi1grid[i][j] = np.ones(ode.I)-np.asarray(ode.final_Xi0())

        fig = plt.figure(figsize=(21,np.ceil(I/3)*10))

        for i in range(I):
            name = self.rxn.phases()[i].name()
            phasein = np.asarray([[name in self.phasenames[i][j] for j in range(len(self.grid.x_range))] for i in range(len(self.grid.y_range))])

            axi = fig.add_subplot(int(np.ceil(I/3)),3,i+1)
            ax = self.setup_axes(axi)
            ax.set_title(name)

            cmap = plt.get_cmap('bwr')
            s = self.scatter(ax,self.grid.xgrid[phasein],self.grid.ygrid[phasein],Xi1grid[phasein][:,i],cmap=cmap)
            fig.colorbar(s,location='left',ax=axi)
            
    @update
    def plot_Ci1(self):
        I = len(self.rxn.phases())
        Ci1grid = np.empty(self.grid.ygrid.shape+(I,))
        for i,P in enumerate(self.grid.y_range):
            for j,x in enumerate(self.grid.x_range):
                ode = self.reconstruct_ode(i,j)
                if ode.sol is not None:
                    Ci1grid[i][j] = np.ones(ode.I)-np.asarray(ode.final_Ci0())

        fig = plt.figure(figsize=(21,np.ceil(I/3)*10))

        for i in range(I):
            name = self.rxn.phases()[i].name()
            phasein = np.asarray([[name in self.phasenames[i][j] for j in range(len(self.grid.x_range))] for i in range(len(self.grid.y_range))])

            axi = fig.add_subplot(int(np.ceil(I/3)),3,i+1)
            ax = self.setup_axes(axi)
            ax.set_title(name)

            cmap = plt.get_cmap('bwr')
            s = self.scatter(ax,self.grid.xgrid[phasein],self.grid.ygrid[phasein],Ci1grid[phasein][:,i],cmap=cmap)
            fig.colorbar(s,location='left',ax=axi)
    
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
        
        def click(event):  
            if event.inaxes==ax:
                
                coords = [(sis[si][i],sjs[si][i]) for si, sc in enumerate(scs) if (contains:=sc.contains(event))[0] for i in contains[1]['ind']]
                if len(coords)>0:
                    ode = self.reconstruct_ode(coords[0][0], coords[0][1])
                    for axis in raxes: 
                        axis.clear()
                        axis.set_visible(True)
                    ode.plot_axes(raxes)
                    ode.apply_labels(raxes)
                    for axis in raxes: axis.set_xlim([0.0,500.])
        
        fig.canvas.mpl_connect("motion_notify_event", hover)
        fig.canvas.mpl_connect("button_press_event", click)



