# -*- coding: UTF-8 -*-

from __future__ import print_function, division

import numpy as np

class CorrFunction(object):
    """
    correlation function data, additional information and manipulation methods

    Parameters
    ----------
    corr : (nvec, nframes) array
        Correlation functions

    std : (nvec, nframes) array
        Correlation function standard deviations

    error : (nvec, nframes) array
        Correlation function standard error of the mean

    info : dict
        Dictionary with information on correlation functions
    """

    def __init__(self, corr=None, std=None, error=None, info=None):
        self.corr        = corr
        self.std         = std
        self.error       = error

        if info is not None:
            self.resid       = info['bondvecinfo']['resid'     ]
            self.resindex    = info['bondvecinfo']['resindex'  ]
            self.resname     = info['bondvecinfo']['resnames'  ]
            self.atomindex   = info['bondvecinfo']['atomindex' ]
            self.atomname    = info['bondvecinfo']['atomnames' ]
            self.element     = info['bondvecinfo']['element'   ]
            self.chain       = info['bondvecinfo']['chain'     ]
            self.bondlength  = info['bondvecinfo']['bondlength']
            self.bondvec     = info['bondvecinfo']['bondvec'   ]
            self.fitgroup    = info['bondvecinfo']['fitgroup'  ]
            try:
                self.fit     = info['bondvecinfo']['fit'       ]
            except KeyError:
                self.fit     = False
            try:
                self.S2direct = np.array(info['bondvecinfo']['S2'])
            except KeyError:
                self.S2direct = None
            self.dt          = info['bondvecinfo']['dt'        ]
            self.topfilename = info['topfilename']
            self.npzfilename = info['npzfilename']
            self.trjfilename = info['trjfilename'] 
            self.frames      = info['frames'     ]
        else:
            self.resid       = None
            self.resindex    = None
            self.resname     = None
            self.atomindex   = None
            self.atomname    = None
            self.element     = None
            self.chain       = None
            self.bondlength  = None
            self.bondvec     = None
            self.fitgroup    = None
            self.fit         = None
            self.dt          = None
            self.topfilename = None
            self.npzfilename = None
            self.trjfilename = None
            self.frames      = None
 
