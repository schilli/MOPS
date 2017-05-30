#!/usr/bin/env python
# This is a demonstration of how to compute bond vector correlation functions.
# A 200 ns Trp-cage (PDB ID: 1L2Y) MD simulation is used as a test case.

from __future__ import print_function, division

import os
import MDorder as mdo

# path to the pdb file and trajectory shipped with MDorder
testpath = '/'.join(mdo.__file__.split('/')[:-1]) + '/test/'

# the correlation functions are stored in a subfolder of the current working directory
savepath = "./MDorder_test_corr_fit"
if not os.path.isdir(savepath):
    os.mkdir(savepath)

# topology and trajectory filenames.
# A list of multiple trajectory filenames can be specified as well
topfilename  =  testpath + "trp-cage.pdb"
trjfilenames = [testpath + "trp-cage.xtc"]

# these are the default selection strings for N and H atoms
Nsel = "name N and not resname PRO and resid >= 1"
Hsel = "name H and not resname PRO and resid >= 1"

# compute bond order correlation functions for a subrajectory length of 10000 ps = 10 ns
mdo.bondvec_corr_batch(topfilename, trjfilenames, savepath, subtrjlength=10000, bondvec=[Nsel, Hsel], fitgroup="backbone")
