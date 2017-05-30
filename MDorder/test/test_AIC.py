#!/usr/bin/env python
# This is a demonstration of how to compute S2 order parameters from bond vector correlation functions.

from __future__ import print_function, division

import sys, os, glob
import MDorder as mdo
import matplotlib.pyplot as plt
 
# the correlation functions are stored in a subfolder of the current working directory
# after running test_corr.py
corrpath = "./MDorder_test_corr"

if not os.path.isdir(corrpath):
    print("No correlation functions found.")
    print("Please run test_corr.py first.")
    sys.exit(1)

# load correlation functions
corrFilenames = glob.glob(corrpath + '/*.zip')
op = mdo.OrderParameter(corrfilenames=corrFilenames)

# predict order parameters
# The best multiexponential model for each correlation function (averaged over subtrajectories) is selected,
# up to 'maxdecays' exponentials
op.estimate_generalLS_modelSelection("generalLSselection", maxdecays=3, weighted=True, nfits=10)

# extract information
S2_all     = op.S2        # array containing order parameters for each exponential decay process. shape = (nresidues, maxdecays)
S2_all_std = op.S2std     # std. dev. corresponding to each op.S2 value (same shape)
S2         = op.S2.max(1) # the max is the traditional S2 order parameter
S2_std     = op.S2std[range(op.S2.shape[0]),op.S2.argmax(1)] # std. dev. corresponding to S2
tau_all    = op.tau       # contains all decay time constants corresponding to the order parameters
tau        = op.tau[range(op.S2.shape[0]),op.S2.argmax(1)] # decay time constants corresponding to S2
tau_std    = op.tau[range(op.S2.shape[0]),op.S2.argmax(1)] # std. dev. of decay time constants corresponding to S2

avgcorr    = op.avgcorr   # correlation function object with averaged correlation functions over all subtrajectories
corr       = avgcorr.corr # numerical correlation functions, array of shape = (nresidues, timeframes)
corrlist   = op.corrlist  # list of correlation functions per subtrajectory

resids    = op.avgcorr.resid[0]    # residue ID of the first residue of the bond vector
residx    = op.avgcorr.resid[0]    # residue index (0-based)
resnames  = op.avgcorr.resname[0]  # residue name
atomnames = op.avgcorr.atomname[0] # atom name

plt.bar(resids, S2, yerr=S2_std)
plt.ylim(0,1)
plt.xlabel('Reisdue Number')
plt.ylabel(r'S$^2$')
plt.show()
