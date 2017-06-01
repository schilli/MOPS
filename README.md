# MDorder
Compute Amide bond S2 order parameters from MD trajectories.

A publication describing the methods used in this program in detail is under way.

# How to Install
MDorder is compatible with python 2 and 3 and runs on linux, maybe on MacOS as well (untested).  
Either simply download and install this package by running the setup script, or:

1. Download and install the [Anaconda python distribution](https://www.continuum.io/downloads "Continuum Analytics Anaconda download")
2. Create a new environment with the [MDTraj package](https://github.com/mdtraj/mdtraj "MDTraj") package installed:  
`conda create --name MDorder mdtraj matplotlib ipython`
3. Activate environment:  
`source activate MDorder`
4. Download and install MDorder:  
`git clone https://github.com/schilli/MDorder.git`  
`cd MDorder`  
`python setup.py intall`
5. To update to a newer version go to the MDorder directory cloned with git and:  
`git pull`  
`python setup.py install`


# Usage

As the computation of the bond vector correlation functions and the computation of order parameters from the correlation functions can both take some time,
their computations have been split into separate steps in the workflow.

You should provide trajectory data that contains only atoms of a single protein, without solvent.

## Standalone Executable

The module automatically installs an executable `MDorder` in the users `$PATH` that will print detailed usage instructions when called with the `-h` flag.

## From a python script
After installation, there will be several files for testing and documentation of the API in your `$PATH`.
These can also be found in the `MDorder/test` subdirectory.
The `which` command will help you to locate them.
They are well documented and demonstrate how to use the API of MDorder from the command line and in your own Python scripts.
In addition, each function and class comes with a docstring describing its usage and the meaning of arguments and return values.
It is a good idea to import the MDorder module in an interactive Ipython session and to read the docstrings with the `?` operator, e.g.:  
```python
In[1]: import MDorder as mdo
In[2]: mdo.OrderParameter?
```

The testing / API documentation scritps available are:
* `test_corr_fit.py`: Computes the internal correlation functions with prior removal of the global rotation of the protein by superposition of backbone atoms.
* `test_corr_nofit.py`: Computes the correlation functions without removal of the global rotation of the protein.
* `test_AIC.py`: Computes order parameters with the general Lipari-Szabo model and selects the best model (i.e., number of exponentials fitted) based on the Aikaike Information Criterion (AIC).
* `test_direct.py` : Computes order parameters with the direct method described in Trbovic et al. Proteins (2008). doi:10.1002/prot.21750
* `test_mean.py` : Computes order parameters as the mean of the internal bond vector correlation function convergence value






