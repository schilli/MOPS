# MDorder
Compute Amide bond S2 order parameters from MD trajectories

## How to Install
Either simply download and install this package by running the setup script, or:

1. Download and install the [Anaconda python distribution](https://www.continuum.io/downloads "Continuum Analytics Anaconda download")
2. Create a new environment with the [MDTraj package](https://github.com/mdtraj/mdtraj "MDTraj") package installed:  
`conda create --name MDorder mdtraj`
3. Activate environment:  
`source activate MDorder`
4. Download and install MDorder:  
`git clone https://github.com/schilli/MDorder.git`  
`cd MDorder`  
`python setup.py intall`
5. To update to a newer version go to the MDorder directory cloned with git and:  
`git pull`  
`python setup.py install`

