"""
Compute amide bond S2 order parameters from MD trajectories. 

The computation is done in two steps:

1. Compute bond vector correlation functions
    mdo.bondvec_corr_batch(topfilename, trjfilenames)

2. Compute S2 order parameters from correlation functions
    op = mdo.OrderParameter(corrfilenames=corrFilenames)
    op.estimate()

Refer to the dosctrings of the individual modules, classes and functions for details.
"""

try:
    # python 3 ?
    from MDorder.MDorder import save_corr, load_corr, bondvec_corr, bondvec_corr_batch
    from MDorder.OrderParameter import OrderParameter
    #from MDorder.CorrFunction import CorrFunction
except ImportError:
    # python 2 ?
    from MDorder import save_corr, load_corr, bondvec_corr, bondvec_corr_batch
    from OrderParameter import OrderParameter 
    #from CorrFunction import CorrFunction

__all__ = []
