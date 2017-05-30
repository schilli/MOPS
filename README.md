# MDorder
Compute Amide bond S2 order parameters from MD trajectories.

A publication describing the methods used in this program in detail is under way.

## How to Install
Either simply download and install this package by running the setup script, or:
MDorder is compatible with python 2 and 3.

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


## Usage

As the computation of the bond vector correlation functions and the computation of order parameters from the correlation functions can both take some time,
their computations have been split into separate steps in the workflow.

# From a python script
After installation, there will be several files for testing in your `$PATH`.
The `which` command will help you to locate them.
They are well documented and demonstrate how to use MDorder from the command line.
In addition, each function and class comes with a docstring describing its use and the meaning of arguments and return values.

* `test_corr_nofit.py`: Computes the correlation functions without removal of the global rotation of the protein.
* `test_AIC.py`: Computes order parameters with the general Lipari-Szabo model and selects the best model (i.e., number of exponentials fitted) based on the Aikaike Information Criterion (AIC).




