# -*- coding: UTF-8 -*-

from __future__ import print_function, division

import sys, os, time, copy
import numpy as np
import matplotlib.pyplot as plt
import LS

from corrFunction import corrFunction

# ============================================================================ #

class OrderParameter(object):
    """
    Bond vector order parameters.

    Parameters
    ----------
    corrfilenames : list
        List with *.zip filenames containing correlation functions.

    converged : boolean, optional, default: True
        Use only converged correlation functions for averaging

    verbose : boolean, optional
        Report progress and other verbose output.

    **kwargs : optional keyword arguments
        All remaining keyword arguments are passed to the order parameter method.
    """

# ==================================== #

    def __init__(self, corrfilenames=None, converged=True, verbose=True, sort=True, **kwargs):
        self.corrfilenames = corrfilenames
        self.converged     = converged
        self.verbose       = verbose

        self.corrlist      = []

        self.method_kwargs = kwargs

        # mean values
        self.avgcorr = None

        # for plotting
        self.figure  = None
        self.axs     = None
        self.cids    = {}   # connection ids for plot event handlers
        self.lines   = []   # line object handels
        self.corrset = 0    # which set of correlation functions from the corrlist should be plottet
        self.corridx = 0    # which correlatin function to plot next
        self.leftButtonPressed = False
        self.plotfit = None

        # load data
        if self.corrfilenames is not None:
            self.load(sort=sort)
 
# ==================================== #

    def load(self, corrfilenames=None, sort=True):
        """
        Load sets of correlation functions.

        Parameters
        ----------
        corrfilenames : list
            List with *.zip filenames containing correlation functions.
        sort: boolean
            if True, sort corrfilenames before loading
 
        """

        if corrfilenames is not None:
            self.corrfilenames = corrfilenames

        # The order in which correlationfunctions are loaded can affect the computation
        # of the mean and standard deviations later on.
        # To ensure reproducability, sort them here!
        if sort:
            self.corrfilenames.sort()

        starttime = time.time()

        # load correlation functions
        self.corrlist = []
        print("Loading {} set{} of correlation functions:".format(len(self.corrfilenames), *['s' if i>1 else '' for i in [len(self.corrfilenames)]]))
        for nf, filename in enumerate(self.corrfilenames):
            corr, corrstd, corrstdmean, info = load_corr(filename)
            self.corrlist.append(corrFunction(corr, corrstd, corrstdmean, info))
            if self.verbose:
                print("\rProgress: {:3.0f}%".format(100.0*nf/len(self.corrfilenames)), end="")
                sys.stdout.flush()

        # compute mean
        self.average_corr()

        # report runtime
        if self.verbose:
            print("\rLoading took: {:.2f} sec.".format(time.time() - starttime))

# ==================================== #

    def average_corr(self, offset=0):
        """
        Compute the average correlation function and its standard deviation

        Parameters:
        -----------
        offset: int
            determines the order in which the correlation functions are averaged
        """
        nvec, nframes = self.corrlist[0].corr.shape
        ncorr         = len(self.corrlist)
        allcorr = np.zeros([ncorr, nvec, nframes], dtype=self.corrlist[0].corr.dtype)

        for i, c in enumerate(range(offset, ncorr+offset)):
            allcorr[i,:,:] = self.corrlist[c % ncorr].corr
        #for c, corr in enumerate(self.corrlist):
        #    allcorr[c,:,:] = corr.corr

        self.avgcorr = corrFunction(corr=allcorr.mean(0), std=allcorr.std(0), error=allcorr.std(0)/allcorr.shape[0]**0.5)
        self.avgcorr.resid       = copy.copy(self.corrlist[0].resid      )
        self.avgcorr.resindex    = copy.copy(self.corrlist[0].resindex   )
        self.avgcorr.resname     = copy.copy(self.corrlist[0].resname    )
        self.avgcorr.atomindex   = copy.copy(self.corrlist[0].atomindex  )
        self.avgcorr.atomname    = copy.copy(self.corrlist[0].atomname   )
        self.avgcorr.element     = copy.copy(self.corrlist[0].element    )
        self.avgcorr.chain       = copy.copy(self.corrlist[0].chain      )
        self.avgcorr.bondlength  = copy.copy(self.corrlist[0].bondlength )
        self.avgcorr.bondvec     = copy.copy(self.corrlist[0].bondvec    )
        self.avgcorr.fitgroup    = copy.copy(self.corrlist[0].fitgroup   )
        self.avgcorr.fit         = copy.copy(self.corrlist[0].fit        )
        self.avgcorr.dt          = copy.copy(self.corrlist[0].dt         )
        self.avgcorr.topfilename = copy.copy(self.corrlist[0].topfilename)
        self.avgcorr.npzfilename = None
        self.avgcorr.trjfilename = None
        self.avgcorr.frames      = None

        return allcorr

# ==================================== #

    def estimate(self, method="mean", converged=True, **kwargs):
        """
        Estimate bond vector order parameters from correlation functions.

        Parameters
        ----------
        method: string, optional
            The method to use for order parameter computation.
            Options are:
                "direct"        If provided, use S2 values directly approximated from the bond vectors as described in: 
                                Trbovic et al. Proteins (2008). doi:10.1002/prot.21750
                "mean"          Use the mean of the final quarter as order parameter
                "single exp"    Fit correlation functions to single exponential
                "double exp"    Fit correlation functions to double exponential
                "LS"            Use the Lipari Szabo model to estimate S2 values
                "extLS"         Use the extended Lipari Szabo method (method 3) from:
                                JPCB 2008, 112, 6155-6158, pubs.acs.org/doi/abs/10.1021/jp077018h

        converged : boolean, optional, default: True
            Use only converged correlation functions for averaging
        """
 
        self.method    = method
        self.converged = converged

        # select order parameter estimation method
        if self.method == "direct":
            self.estimate_direct(**kwargs)
        elif self.method == "mean":
            self.estimate_mean(converged=self.converged, **kwargs)
        elif self.method == "generalLS":
            self.estimate_generalLS(**kwargs)
        elif self.method == "generalLSselection":
            self.estimate_generalLS_modelSelection(**kwargs)
        else:
            print("Order parameter estimation method unknown: {}".format(self.method))
            sys.exit(1)
 
# ==================================== #

    def estimate_direct(self):
        """
        Compute mean bond vector order parameters from direct estimates as described in:
        Trbovic et al. Proteins (2008). doi:10.1002/prot.21750
        """
        if self.corrlist[0].S2direct is not None:
            self.S2all = np.zeros([self.corrlist[0].corr.shape[0], len(self.corrlist)], dtype=np.float)
            for n, corrfun in enumerate(self.corrlist):
                self.S2all[:,n] = corrfun.S2direct
            self.S2mean  = self.S2all.mean(1)
            self.S2std   = self.S2all.std(1)
            self.S2error = self.S2all.std(1) / self.S2all.shape[1]**0.5 
        else:
            print("Direct estimate of S2 not possible. Data not present in correlation functions.")

# ==================================== #

    def estimate_mean(self, converged=True, diffThreshold=0.02, stdThreshold=0.02):
        """
        Estimate bond vector order parameters as the mean of the last quarter of the correlation function.

        Parameters
        ----------
        converged : boolean, optional
            If True, only use converged correlation functions.
        diffThreshold : float, optional
            If the quarters 3 and 4 of the correlation function differ more than this threshold,
            they are not considered converged
        stdThreshold : float, optional
            If the mean standard deviation of quarter 3 and 4 is larger than this threshold,
            they are not considered converged.
        """

        self.S2all         = np.zeros([self.corrlist[0].corr.shape[0], len(self.corrlist)], dtype=np.float)
        self.S2convergence = np.zeros_like(self.S2all, dtype=np.bool)

        # compute S2 values for each correlation function and judge convergence
        for c, corrfun in enumerate(self.corrlist):
            mean, convergence = self.estimate_mean_single(corrfun, diffThreshold=diffThreshold, stdThreshold=stdThreshold)
            self.S2all[:,c]         = mean
            self.S2convergence[:,c] = convergence

        # compute global S2 as mean of individiual S2
        if not converged:
            self.S2mean  = self.S2all.mean(1)
            self.S2std   = self.S2all.std(1)
            self.S2error = self.S2all.std(1) / self.S2all.shape[1]**0.5
        else:
            self.S2mean  = np.zeros(self.S2all.shape[0])
            self.S2std   = np.zeros(self.S2all.shape[0])
            self.S2error = np.zeros(self.S2all.shape[0])
            for n in range(self.S2all.shape[0]):
                self.S2mean[n]  = self.S2all[n,self.S2convergence[n,:]].mean()
                self.S2std[n]   = self.S2all[n,self.S2convergence[n,:]].std()
                self.S2error[n] = self.S2std[n] / self.S2convergence[n,:].sum()**0.5

        self.S2nconverged = self.S2convergence.sum(1)

        # compute global S2 from average correlation functions
        self.S2avg, self.S2avgConvergence = self.estimate_mean_single(self.avgcorr, diffThreshold=diffThreshold, stdThreshold=stdThreshold)

# ==================================== #

    def estimate_mean_single(self, corrfun, converged=True, diffThreshold=0.02, stdThreshold=0.02):
        """
        Estimate bond vector order parameters for a single set of correlation functions.
        
        Parameters
        ----------
        corrfun : (nvec, nframes) array
            Single set of correlation functions
        other parameters:
            see function "estimate_mean()"
        """

        length            = corrfun.corr.shape[1]
        quarter           = length / 4
        thirdQuarter      = corrfun.corr[:,2*quarter:3*quarter]
        fourthQuarter     = corrfun.corr[:,3*quarter:4*quarter]
        fourthQuarterMean = fourthQuarter.mean(1)

        difference        = abs(thirdQuarter.mean(1) - fourthQuarterMean)
        stdev             = (thirdQuarter.std(1) + fourthQuarter.std(1)) / 2
        convergence       = np.logical_and(difference < diffThreshold, stdev < stdThreshold)

        return fourthQuarterMean, convergence 

# ==================================== #

    def estimate_generalLS(self, n, fast=True, internal=False, weighted=False, **kwargs):

        dt      = self.avgcorr.dt
        ncorr   = self.avgcorr.corr.shape[0]
        nframes = self.avgcorr.corr.shape[1]
        t       = np.linspace(0, dt*nframes, nframes)
        firstf  = 0
        if fast:
            firstf = 1

        self.para = []
        self.S2   = np.zeros([ncorr, n])
        self.tau  = np.zeros([ncorr, n])

#        nS   = n
#        ntau = n 
#        startvals  = []
#        startvals.append(2*n*[1.0])
#        startvals.append(list(np.linspace(0.1, 0.9, nS)) + ntau*[1.0])
#        startvals.append(nS*[1.0] + list(np.logspace(0, ntau-1, ntau))[-1::-1])
#        startvals.append(list(np.linspace(0.1, 0.9, nS)) + list(np.logspace(0, ntau-1, ntau))[-1::-1])
#        self.lsq = np.zeros([ncorr, len(startvals)])

        for nc in range(ncorr):
            if weighted:
                self.ls = LS.LS(t[firstf:], self.avgcorr.corr[nc,firstf:], sigma=self.avgcorr.std[nc,firstf:])
            else:
                self.ls = LS.LS(t[firstf:], self.avgcorr.corr[nc,firstf:])
            p = self.ls.fit(n, fast=fast, internal=internal, **kwargs)

#            for i, p0 in enumerate(startvals):
#                p = self.ls.fit(n, fast=fast, internal=internal, p0=p0, **kwargs)
#                self.lsq[nc,i] = p["lsq"]
#
#            print("{:3d} done".format(nc))

            self.para.append(p)
            self.S2[nc,:]  = p["S"]
            self.tau[nc,:] = p["tau"]

# ==================================== #

    def check_overfitting(self, parameters, mintauratio=2.0, minS2diff=0.01):
        """
        Check generalLS parameters for overfitting
        """
        overfitted = False

        if not parameters['success']:
            return True

        if len(parameters['tau']) > 1:
            minratio = min(parameters['tau'][:-1] / parameters['tau'][1:])
            if minratio < mintauratio:
                overfitted = True
#                print("tau overfitted:", parameters['tau'])

        if len(parameters['S']) > 1:
            mindiff = min(abs(parameters['S'][:-1] - parameters['S'][1:]))
            if mindiff < minS2diff:
                overfitted = True
#                print("S overfitted:", parameters['S'])

        if np.isnan(parameters['p']).sum() > 0:
            overfitted = True

        return overfitted


# ==================================== #

    def estimate_generalLS_modelSelection(self, fast=True, internal=False, weighted=False, maxdecays=int(1e3), nfits=1, ncorr=None, **kwargs):
        """
        nfits
        The number of fits to perform per amino acid.
        The order of the correlation function for the mean correlation function is randomized for each fit.
        As fitting is an ill-posed problem, we get better estimates in this way
        """

        randomstate = np.random.get_state()
        np.random.seed(23)

        dt      = self.avgcorr.dt
        if ncorr is None:
            ncorr   = self.avgcorr.corr.shape[0]
#        ncorr   = 10
        nframes = self.avgcorr.corr.shape[1]
        t       = np.linspace(0, dt*nframes, nframes)
        firstf  = 0
        if fast:
            firstf = 1

        self._AIC      = np.zeros([ncorr, nfits, maxdecays])
        self._lsq      = np.zeros([ncorr, nfits, maxdecays])
        self._para     = np.zeros([ncorr, nfits, maxdecays, 2*maxdecays])
        self._S2       = np.zeros([ncorr, nfits, maxdecays,   maxdecays])
        self._tau      = np.zeros([ncorr, nfits, maxdecays,   maxdecays])
        self._success  = np.zeros([ncorr, nfits, maxdecays], dtype=np.bool)
        self._paralist = [[[] for nfit in range(nfits)] for nc in range(ncorr)]

        self.para = []

        def make_progress_msg(progress_percent, ETA):
            progress_msg = "Progress: {:3.0f}% ".format(progress_percent)
            if ETA > 60*60*24:
                ETA_msg = "(ETA: {:5.1f} days)".format(ETA/(60*60*24))
            elif ETA > 60*60:
                ETA_msg = "(ETA: {:5.1f} hours)".format(ETA/(60*60))
            elif ETA > 60:
                ETA_msg = "(ETA: {:5.1f} minutes)".format(ETA/(60))
            else:
                ETA_msg = "(ETA: {:5.0f} seconds)".format(ETA)
            return progress_msg + ETA_msg


        ETA = 0.0
        progress_msg = make_progress_msg(0.0, ETA)
        print(progress_msg, end="")
        sys.stdout.flush()

        starttime = time.time()
        for nc in range(ncorr):

            for nfit in range(nfits):
                print(len(progress_msg)*'\b' + len(progress_msg)*' ' + len(progress_msg)*'\b', end="")
                progress_percent = 100.0*(nc*nfits+nfit)/(nfits*ncorr)
                runtime          = time.time() - starttime
                if progress_percent > 0:
                    ETA          = runtime * (100-progress_percent) / progress_percent
                progress_msg = make_progress_msg(progress_percent, ETA)
                print(progress_msg, end="")
                sys.stdout.flush()

                # compute new average correlation function
                np.random.shuffle(self.corrlist)
                self.average_corr()

                # set up Lipari Szabo fitter
                if weighted:
                    self.ls = LS.LS(t[firstf:], self.avgcorr.corr[nc,firstf:], sigma=self.avgcorr.std[nc,firstf:])
                else:
                    self.ls = LS.LS(t[firstf:], self.avgcorr.corr[nc,firstf:]) 

                # fit for all correlation functions, correlation function shuffles and number of decays
                for ndecays in range(1,maxdecays+1):
                    decay_ndx = ndecays - 1
                    p = self.ls.fit(ndecays, fast=fast, internal=internal, **kwargs)
                    self._paralist[nc][nfit].append(p)
                    self._AIC    [nc, nfit, decay_ndx]             = p['AIC']
                    self._lsq    [nc, nfit, decay_ndx]             = p['lsq']
                    self._para   [nc, nfit, decay_ndx, :2*ndecays] = p['p']
                    self._S2     [nc, nfit, decay_ndx,   :ndecays] = p['S']
                    self._tau    [nc, nfit, decay_ndx,   :ndecays] = p['tau']
                    self._success[nc, nfit, decay_ndx]             = p['success']

        print(len(progress_msg)*'\b' + len(progress_msg)*' ' + len(progress_msg)*'\b', end="")
        progress_msg = make_progress_msg(100.0, 0.0)
        print(progress_msg)

        # reset random number generator state
        np.random.set_state(randomstate)

        # compute probabilities for each model to be the best
        self._probability  = np.exp((self._AIC.min(2, keepdims=True) - self._AIC)/2)
        self._probability /= self._probability.sum(2, keepdims=True)
        self._meanprob     = self._probability.mean(1)

        # select best model based on AIC
        self._bestmodel = self._meanprob.argmax(1)
        self._bestprob  = self._meanprob.max(1)
        self._goodfits  = self._probability.argmax(2) - self._bestmodel.reshape([self._bestmodel.shape[0],1]) == 0

        # store S2 and tau of the best model for each residue
        self.S2      = np.array([self._S2 [nc,self._goodfits[nc,:],self._bestmodel[nc],:].mean(0) for nc in range(ncorr)])
        self.S2std   = np.array([self._S2 [nc,self._goodfits[nc,:],self._bestmodel[nc],:].std (0) for nc in range(ncorr)])
        self.tau     = np.array([self._tau[nc,self._goodfits[nc,:],self._bestmodel[nc],:].mean(0) for nc in range(ncorr)])
        self.taustd  = np.array([self._tau[nc,self._goodfits[nc,:],self._bestmodel[nc],:].std (0) for nc in range(ncorr)]) 
        self.ndecays = self._bestmodel + 1
        self.prob    = self._bestprob

        # store list of one original fitting result for each residue
        self.para = []
        for nc in range(ncorr):
#            p = {}
#            p['S']   = np.array([1, 0.9, 0.8])
#            p['tau'] = np.array([3000, 1000, 50])
#            p['p']   = np.array([1, 0.9, 0.8, 3000, 1000, 50])
#            self.para.append(p)
            clearestfit = self._probability[nc,:,self._bestmodel[nc]].argmax()
            self.para.append(self._paralist[nc][clearestfit][self._bestmodel[nc]])
            self.para[-1]['S'  ] = self.S2 [nc,:self.ndecays[nc]]
            self.para[-1]['tau'] = self.tau[nc,:self.ndecays[nc]]
            self.para[-1]['p'  ] = np.concatenate((self.para[-1]['S'], self.para[-1]['tau']))


# ==================================== #

    def plot_all(self, savepath, internal=False):

        fig = plt.figure()
        axs = fig.add_subplot(111)

        nc, npoints  = self.avgcorr.corr.shape
        dt      = self.avgcorr.dt
        t       = np.linspace(0, dt*npoints, npoints)

        if not os.path.isdir(savepath):
            os.mkdir(savepath)

        for c in range(nc):
            corr    = self.avgcorr.corr[c,:]
            corrstd = self.avgcorr.std[c,:]
            correrr = self.avgcorr.error[c,:]
            para    = self.para[c]
            ls      = LS.LS(t, corr)
            ls.internal = internal
            corrfit = ls.generalLS(para['p'])

            coeffs  = ",  ".join([r"$S_{{{}}}^{{2}}={:.2f}$".format(i,S) for i,S in enumerate(para['S'])])
            coeffs += ",  "
            coeffs += ",  ".join([r"$\tau_{{{}}}={:.1f}$".format(i,tau) for i,tau in enumerate(para['tau'])]) 

            lines = []
            lines += axs.plot(t, corr   , 'b', label="Mean correlation function over {} samples".format(len(self.corrlist)), lw=2)
            lines += axs.plot(t, corrfit, 'g', label="fit: {}".format(coeffs), lw=2)
            #lines += axs.fill_between(t, corr+1.96*correrr, corr-1.96*correrr, alpha=0.4, color='r', label="95% Confidence Interval")

            resid   = self.avgcorr.resid[0][c]
            resname = self.avgcorr.resname[0][c]
            atom1   = self.avgcorr.atomname[0][c]
            atom2   = self.avgcorr.atomname[1][c]

            axs.set_ylim(min(0, corr.min()), max(1, corr.max()))
            axs.set_title("{} {} ({}-{})".format(resname, resid, atom1, atom2))
            axs.set_ylabel("correlation")
            axs.set_xlabel("time [ps]")
            lgd = axs.legend(loc="upper right", bbox_to_anchor=(1.2,1)) 

            filename = savepath + "/corr_{:03d}.jpg".format(resid)
            print(filename)
            plt.savefig(filename, dpi=100, papertype='a4', orientation='landscape', bbox_inches='tight', bbox_extra_artists=(lgd,))

            while len(lines) > 0:
                line = lines.pop()
                line.remove()  



        plt.close()

# ==================================== #

    def plot_corr(self, event=None, corrset=None, fit=None):
        """
        Plot correlation functions.
        
        Parameters
        ----------
        event : matplotlib mouse event
            Only there to update the plot on various matplotlib events, like mouse button presses
        corrset: integer, optional
            Which set of correlation functions to use for plotting.
            If None, plot the mean.
        fit: string
            What fit should be plotted? Options are:
                "double exp"
                "LS"
        """

        if fit is not None:
            self.plotfit = fit

        if self.figure is None:
            # create figure and axis
            self.figure = plt.figure()
            self.axs    = self.figure.add_subplot(111)  
            self.cids['button_press'  ] = self.figure.canvas.mpl_connect('button_press_event',   self._onclick)
            self.cids['button_release'] = self.figure.canvas.mpl_connect('button_release_event', self._onrelease)
            self.cids['motion'        ] = self.figure.canvas.mpl_connect('motion_notify_event',  self._onmove)
            self.cids['scroll'        ] = self.figure.canvas.mpl_connect('scroll_event',         self._onscroll)
            self.cids['close'         ] = self.figure.canvas.mpl_connect('close_event',          self._onclose)

            if corrset is None:
                self.corrset  = 0
                self.plotMean = True
            else:
                self.corrset  = corrset
                self.plotMean = False
            self.corridx = 0
            self.lines   = []

            if self.avgcorr is None:
                self.average_corr()

        # remove old data
        while len(self.lines) > 0:
            line = self.lines.pop()
            line.remove() 

        start_idx = 0

        # plot data
        if self.plotMean:
            corrFun     = self.avgcorr
            xdata       = np.linspace(0, corrFun.corr.shape[1] * corrFun.dt, corrFun.corr.shape[1])    
            self.lines += self.axs.plot(xdata[start_idx:], self.avgcorr.corr[self.corridx,start_idx:], 'b', label="Mean correlation function over {} samples".format(len(self.corrlist)), lw=2)
            #self.lines.append(self.axs.fill_between(xdata[start_idx:], self.avgcorr.corr[self.corridx,start_idx:]+self.avgcorr.std[self.corridx,:],
            #                                               self.avgcorr.corr[self.corridx,start_idx:]-self.avgcorr.std[self.corridx,:],
            #                                               alpha=0.4, color='r')) 
            self.lines.append(self.axs.fill_between(xdata[start_idx:], self.avgcorr.corr[self.corridx,start_idx:]+1.96*self.avgcorr.error[self.corridx,start_idx:],
                                                           self.avgcorr.corr[self.corridx,start_idx:]-1.96*self.avgcorr.error[self.corridx,start_idx:],
                                                           alpha=0.4, color='r', label="95% Confidence Interval")) 

            try:
                #coeffs = ",  ".join(["{:.1f}".format(n) for n in self.poptavg[self.corridx,:]])
                coeffs  = ",  ".join([r"$S_{{{}}}^{{2}}={:.2f}$".format(i,S) for i,S in enumerate(self.para[self.corridx]['S'])])
                coeffs += ",  "
                coeffs += ",  ".join([r"$\tau_{{{}}}={:.1f}$".format(i,S) for i,S in enumerate(self.para[self.corridx]['tau'])])
            except AttributeError as e:
                print("Did you forget to fit the data?")
                raise e
            xdata  = np.linspace(0, corrFun.corr.shape[1] * corrFun.dt, 10*corrFun.corr.shape[1])  
            if self.plotfit == "generalLS":
                self.lines += self.axs.plot(self.ls.t, self.ls.generalLS(self.para[self.corridx]['p']), label="fit: {}".format(coeffs), color='g', lw=2)
            elif self.plotfit in ["single exp", "double exp"]:
                self.lines += self.axs.plot(xdata[start_idx:], n_exp(xdata[start_idx:], *self.poptavg[self.corridx,:]), label="fit: {}".format(coeffs), color='g', lw=2)
            elif self.plotfit == "LS":
                self.lines += self.axs.plot(xdata[start_idx:], LipariSzabo(xdata[start_idx:], *self.poptavg[self.corridx,:]), label="fit: {}".format(coeffs), color='g', lw=2)
            elif self.plotfit == "extLS":
                self.lines += self.axs.plot(xdata[start_idx:], extLipariSzabo(xdata[start_idx:], *self.poptavg[self.corridx,:]), label="fit: {}".format(coeffs), color='g', lw=2)
            elif self.plotfit == "CloreLS":
                coeffs = ",  ".join(["{:.1f}".format(n) for n in self.poptavg[self.corridx,:]])
                self.lines += self.axs.plot(xdata[start_idx:], CloreLS(xdata[start_idx:], *self.poptavg[self.corridx,:]), label="fit: {}".format(coeffs), color='g', lw=2) 
            elif self.plotfit == "extCloreLS":
                label = r"fit: S$^2$={:.2f},  S$_f$$^2$={:.2f},  S$_i$$^2$={:.2f},  $\tau_f$={:.1f},  $\tau_e$={:.1f},  $\tau_m$={:.1f}".format(*self.poptavg[self.corridx,:])
                self.lines += self.axs.plot(xdata[start_idx:], extCloreLS(xdata[start_idx:], *self.poptavg[self.corridx,:]), label=label, color='g', lw=2)  

        else:
            corrFun     = self.corrlist[self.corrset]
            xdata       = np.linspace(0, corrFun.corr.shape[1] * corrFun.dt, corrFun.corr.shape[1])   
            self.lines += self.axs.plot(xdata[start_idx:], corrFun.corr[self.corridx,start_idx:], 'b', lw=2)
            try:
                coeffs = " ".join(["{:.1e}".format(n) for n in self.poptavg[self.corridx,:]])
            except AttributeError as e:
                print("Did you forget to fit the data?")
                raise e 
            xdata  = np.linspace(0, corrFun.corr.shape[1] * corrFun.dt, 10*corrFun.corr.shape[1])  
            if self.plotfit in ["single exp", "double exp"]:
                self.lines += self.axs.plot(xdata[start_idx:], n_exp(xdata[start_idx:], *self.poptavg[self.corridx,start_idx:, self.corrset]), label="fit: {}".format(coeffs), color='g', lw=2)
            elif self.plotfit == "LS":
                self.lines += self.axs.plot(xdata[start_idx:], LipariSzabo(xdata[start_idx:], *self.poptavg[self.corridx,:, self.corrset]), label="fit: {}".format(coeffs), color='g', lw=2)
            elif self.plotfit == "extLS":
                self.lines += self.axs.plot(xdata[start_idx:], extLipariSzabo(xdata[start_idx:], *self.poptavg[self.corridx,:]), label="fit: {}".format(coeffs), color='g', lw=2) 

        # set axis limits
        self.axs.set_ylim(min(0, corrFun.corr.min()), max(1, corrFun.corr.max()))

        # plot scrollbar
        xmin, xmax = self.axs.get_xlim()
        self.lines += self.axs.plot([xmin, xmax], 2*[0.98], 'k', linewidth=2)
        self.lines += self.axs.plot(xmin + self.corridx*(xmax-xmin)/(corrFun.corr.shape[0]-1), 0.98, 'sk', markersize=15)

        # annotate plot
        self.axs.set_title("{} {} ({}-{})".format(corrFun.resname[0][self.corridx], corrFun.resid[0][self.corridx], corrFun.atomname[0][self.corridx], corrFun.atomname[1][self.corridx]))
        self.axs.set_ylabel("correlation")
        self.axs.set_xlabel("time [ps]")
        self.axs.legend(loc="upper right")
        #self.axs.legend(loc="lower left")

        self.figure.canvas.draw()
        #plt.savefig("plots/extCloreLS/corr_{:03d}.jpg".format(corrFun.resid[0][self.corridx]), dpi=100, papertype='a5', orientation='landscape', bbox_inches='tight')
        plt.show()
 
# ==================================== #

    def _onclick(self, event):
        if event.button == 1:
            self.leftButtonPressed = True
            try:
                xmin, xmax = self.axs.get_xlim()
                try:
                    ncorr = self.corrlist[self.corrset].corr.shape[0]
                except IndexError:
                    ncorr = self.avgcorr.corr.shape[0] 
                self.corridx = int(np.round((ncorr-1) * (event.xdata - xmin) / (xmax - xmin)))
                if self.corridx < 0:
                    self.corridx = 0
                elif self.corridx >= ncorr:
                    self.corridx = ncorr - 1
                self.plot_corr()
            except:
                pass 

# ==================================== #

    def _onrelease(self, event):
        if event.button == 1:
            self.leftButtonPressed = False 

# ==================================== #

    def _onmove(self, event):
        if self.leftButtonPressed:
            try:
                xmin, xmax = self.axs.get_xlim()
                try:
                    ncorr = self.corrlist[self.corrset].corr.shape[0]
                except IndexError:
                    ncorr = self.avgcorr.corr.shape[0]
                self.corridx = int(np.round((ncorr-1) * (event.xdata - xmin) / (xmax - xmin)))
                if self.corridx < 0:
                    self.corridx = 0
                elif self.corridx >= ncorr:
                    self.corridx = ncorr - 1 
                self.plot_corr()
            except:
                pass 

# ==================================== #

    def _onscroll(self, event):
        try:
            ncorr = self.corrlist[self.corrset].corr.shape[0]
        except IndexError:
            ncorr = self.avgcorr.corr.shape[0]
        self.corridx -= int(event.step)
        if self.corridx < 0:
            self.corridx = ncorr - 1
        elif self.corridx >= ncorr:
            self.corridx = 0

        self.plot_corr()
 
# ==================================== #
    
    def _onclose(self, event):
        self.figure  = None
        self.axs     = None
        self.cids    = {}
        self.lines   = []
             
