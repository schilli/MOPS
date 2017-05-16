"""
General Lipari Szabo model fits
"""

from __future__ import print_function

import sys
import numpy as np
from scipy.optimize import minimize

# ============================================================================ #

class LS(object):
    """
    Fit data to the general Lipari Szabo model.
    """

    # ================================ #

    def __init__(self, t, C, sigma=None):
        """
        Parameters
        ----------
        t : (N, ) array
            timeseries
        C : (M, N) array
            single correlation function or set of correlation functions
        """

        self.t = t 
        self.C = C
        if sigma is not None:
            self.sigma  = sigma
            self.sigma2 = sigma**2
        else:
            self.sigma  = np.ones_like(C)
            self.sigma2 = self.sigma
        self.internal = False

    # ================================ #

    def fit(self, n, fast=True, internal=False, taum=None, taumax=None, p0=None, ftol=1e-12, eps=1.4901161193847656e-08, maxiter=1000, verbose=False):
        """
        Fit the data and return optimal parameters

        taum sets a fixed global rotation time
        """

        nS   = n
        ntau = n
        self.internal = internal
        if self.internal:
            ntau -= 1

        if p0 is None or len(p0) != nS + ntau:
            # set initial parameters
            p0  = list(np.linspace(0.1, 0.9, nS))
            p0 += ntau*[1.0]
            #p0 += list(np.logspace(0, ntau-1, ntau))[-1::-1]

        # set bounds on parameters
        taumin = 1*(self.t[1] - self.t[0])
        bounds  = nS   * [(0,1)]    # S
        bounds += ntau * [(taumin,None)] # tau


        # set inequality constraints on S parameters
        constraints = []
        for idx in range(nS-1):
            constraints += [{"type": "ineq", "fun":  lambda x,i=idx: x[i+1] - x[i]}]  # Si+1 > Si

        # set inequality constraints on tau parameters
        for idx in range(nS, nS+ntau-1):
            constraints += [{"type": "ineq", "fun":  lambda x,i=idx: x[i] - x[i+1]}]  # tau_{i+1} > tau_{i}

        # if not fast, constrain Sf to 1.0
        if not fast:
            constraints += [{"type": "ineq", "fun":  lambda x,i=nS: x[i-1] - 1.0}] # constrain Sf to 1.0
#        else:
#            constraints += [{"type": "ineq", "fun":  lambda x,i=nS: x[i-1] - self.C[0]}] # constrain Sf to first datapoint

        if not internal and taum is not None:
            constraints += [{"type": "eq", "fun":  lambda x,: x[nS] - taum}] # constrain taum to user provided value
        elif taumax is not None:
            bounds[nS] = (0,taumax)



        # fit
        options = {"ftol": ftol, "eps": eps, "maxiter": maxiter, "disp":verbose}
        res = minimize(self.generalLS_obj, x0=p0, bounds=bounds, constraints=constraints, method="SLSQP", options=options) # L-BFGS-B, TNC, COBYLA, SLSQP

        # extract optimal parameters
        result = {}
        result["p"]          = res.x # all parameters for function input
        result["S"]          = res.x[:nS]
        result["tau"]        = res.x[nS:]
        result["lsq"]        = res.fun
        result["success"]    = res.success
        result["iterations"] = res.nit
        result["status"]     = res.status
        result["message"]    = res.message
        result["AIC"]        = AIC(RSS=result['lsq'], n=self.C.shape[0], k=result["p"].shape[0]+1) # k=len(p) + 1,
                                                                                                   # because noise has a parameter has well 
                                                                                                   # (https://en.wikipedia.org/wiki/Akaike_information_criterion)
        if self.internal:
            result["tau"] = np.array(list(result["tau"]) + [float('inf')])

        return result

    # ================================ #

    def generalLS_internal(self, S, tau):
        """
        Evaluate the general internal Lipari Szabo model.
        """

        C_I = S[0]
        for i, taui in enumerate(tau):
            C_I += (S[i+1] - S[i]) * np.exp(-self.t/taui)
        return C_I

    # ================================ #

    def generalLS(self, p):
        """
        General LS model
        
        Parameters
        ----------
        t : (N, ) array
            timeseries
        p : (2*n, ) array
        """
        n   = len(p)+1
        S   = p[:n/2]
        tau = p[n/2:] 

        if self.internal:
            C_I = self.generalLS_internal(S, tau)
            return C_I
        else:
            C_I = self.generalLS_internal(S, tau[1:])
            C   = C_I * np.exp(-self.t/tau[0])
            return C
        
    # ================================ #

    def generalLS_obj(self, p):
        """
        Objective function that the minimizer is going to minimize.
        Returns the sum of squares of the residual of data and fit.
        """

        C             = self.generalLS(p)
        diff          = C - self.C
        self.squares  = diff**2
        self.weighted = self.squares / self.sigma
        sumsquares    = self.weighted.sum()
        return sumsquares

# ==================================== #

def AIC(RSS, n, k):
    """
    Akaike information criterion according to:
    https://en.wikipedia.org/wiki/Akaike_information_criterion#Comparison_with_least_squares

    Parameters
    ----------
    RSS, double
        Residual sum of squares: sum^{n}_{i=1} (y_i - f(x_i))**2
        with n the number of datapoints and f() the model
    n, integer
        number of datapoints
    k, integer
        number of parameters of the model f()
    """

    return 2*k + n*np.log(RSS)
 
