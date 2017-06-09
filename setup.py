#!/usr/bin/env python

from distutils.core import setup

setup(name='MOPS',
      version='0.1',
      description='Molecular Order Parameters S2 -- Compute Amide bond S2 order parameters from MD trajectories',
      author='Oliver Schillinger',
      author_email='oliver@schillinger-ttp.de',
      url='https://github.com/schilli/MOPS',
      packages=['MOPS'],
      package_data={'MOPS': ['demo/*']},
      scripts=['MOPS/demo/demo_corr_fit.py',
          'MOPS/demo/demo_corr_nofit.py',
          'MOPS/demo/demo_AIC.py',
          'MOPS/demo/demo_mean.py',
          'MOPS/demo/demo_direct.py',
          'MOPS/bin/MOPS']
     )

