"""
Compute amide bond S2 order parameters from MD trajectories. 

The computation is done in two steps:

1. Compute bond vector correlation functions
    mops.bondvec_corr_batch(topfilename, trjfilenames)

2. Compute S2 order parameters from correlation functions
    op = mops.OrderParameter(corrfilenames=corrFilenames)
    op.estimate()

Refer to the dosctrings of the individual modules, classes and functions for details.
"""

__version__ = (1, 0)

try:
    # python 3 ?
    from MOPS.MOPS import save_corr, load_corr, bondvec_corr, bondvec_corr_batch
    from MOPS.OrderParameter import OrderParameter
    #from MOPS.CorrFunction import CorrFunction
except ImportError:
    # python 2 ?
    from MOPS import save_corr, load_corr, bondvec_corr, bondvec_corr_batch
    from OrderParameter import OrderParameter 
    #from CorrFunction import CorrFunction

__all__ = []
