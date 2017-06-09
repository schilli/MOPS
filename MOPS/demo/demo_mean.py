#!/usr/bin/env python
# This is a demonstration of how to compute S2 order parameters from bond vector correlation functions.
# The S2 estimation is done with the method described in:
# Trbovic et al. Proteins (2008). doi:10.1002/prot.21750

from __future__ import print_function, division

import sys, os, glob
import MOPS as mops
import matplotlib.pyplot as plt
 
# the correlation functions are stored in a subfolder of the current working directory
# after running test_corr.py
corrpath = "./MOPS_test_corr_fit"

if not os.path.isdir(corrpath):
    print("No correlation functions found.")
    print("Please run test_corr_fit.py first.")
    sys.exit(1)

# load correlation functions
corrFilenames = glob.glob(corrpath + '/*.zip')
op = mops.OrderParameter(corrfilenames=corrFilenames)

# predict order parameters, take only converged correlation functions into account
op.estimate("mean", converged=True)
 
# extract information
S2         = op.S2mean
S2_std     = op.S2std
S2_err     = op.S2error # = S2.std / <number subtrajectories>

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
