import io
import numpy as np
from contextlib import redirect_stdout, redirect_stderr
import time

from ..base import *
from .base import BasePDReactiveODE

SCIPY_AVAIL = True
try: 
    from scipy.integrate import solve_ivp
except ImportError:
    SCIPY_AVAIL = False

class ScipyPDReactiveODE(BasePDReactiveODE):
    def solve(self,T,p,mi0,Cik0,end,**kwargs):
        super().set_initial_params(T,p,mi0,Cik0,**kwargs)

        method = kwargs.get('method', 'BDF')
        rtol   = kwargs.get('rtol', 1.e-5)
        atol   = kwargs.get('atol', 1.e-9)

        u0 = np.concatenate((mi0,Cik0))

        try:
          so = io.StringIO()
          se = io.StringIO()
          with redirect_stdout(so), redirect_stderr(se):
            tic = time.perf_counter()
            self.sol = solve_ivp(self.rhs, [0,end], u0, dense_output=True, \
                                 method=method, rtol=rtol, atol=atol, jac=self.jac)
            toc = time.perf_counter()
            self.stime = toc-tic
          self.stdout = so.getvalue()
          self.stderr = se.getvalue()
          flag = self.rxn.check_coder_error()
          if flag!=0:
            self.sol = None
            self.excstr = repr(flag)
            self.rxn.reset_coder_error()
        except Exception as e:
          self.sol = None
          self.excstr = repr(e)

if not SCIPY_AVAIL:
    class ScipyPDReactiveODE(BasePDReactiveODE):
        def __init__(self,rxn):
            raise Exception("scipy integrate not available")


