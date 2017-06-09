# -*- coding: UTF-8 -*-
"""
Compute NMR observables from MD trajectories

The usage of MPI is optional.
The MPI parallelized batch function is only defined if mpi4py is installed.
"""

from __future__ import print_function, division

import sys, os, time, string
import mdtraj as md
import numpy as np
import zipfile, bz2
try:
    # python 2 ? 
    import cPickle as pickle
except ImportError:
    # python 3 ? 
    import pickle

try:
    # python 3 ? 
    from MOPS.CorrFunction import CorrFunction
except ImportError:
    # python 2 ? 
    from CorrFunction import CorrFunction
 
try:
    from mpi4py import MPI
    gotMPI = True
except ImportError:
    gotMPI = False


# ============================================================================ #
 

def _fit_trj(trj, fitframe=0, ref=None, fitgroup="name CA", parallel=True):
    """
    Fit MD trajectory to reference based on the specified set of atoms.
    Fitting is done in place.

    Parameters
    ----------
    trj : mdtraj trajectory
        Trajectory to fit.

    fitframe : int, optional
        Frame index of reference frame.

    ref : mdtraj trajectory
        Reference MDTRAJ trajectory.
        If none, the first trajectory is used as a reference.

    fitgroup : string or integer array
        mdtraj selection string or integer array of atom indices
        to select a set of atoms for RMSD fit.

    parallel : boolean, optional
        Toggle parallel fitting algorithm of mdtraj
    """
    if type(ref) == type(None):
        ref = trj

    try:
        fit_atomndx = trj.topology.select(fitgroup)
    except:
        fit_atomndx = fitgroup

    trj.superpose(ref, fitframe, atom_indices=fit_atomndx, parallel=parallel)
 

# ============================================================================ #

def _vec_atoms(trj, bondvec=None):
    """
    Compute bond atom indices.

    Parameters
    ----------
    trj : mdtraj trajectory
        Trajectory for topology information.

    bondvec : 2 element list/tuple, optional
        Both elements specify atoms to compute the bond vector from, either as mdtraj selection
        strings or as arrays with integer atom indices. Defaults to NH protein backbone bond vectors.

    Returns
    -------
    bondvec_ndx : (nvec, 2) array
        Integer array of all atom indices for bond vector atoms
    """

    # default option
    if type(bondvec) == type(None):
        selectiontext1   = "name N and not resname PRO and resid >= 1"
        selectiontext2   = "name H and not resname PRO and resid >= 1"
        atom1_ndx        = trj.topology.select(selectiontext1)
        atom2_ndx        = trj.topology.select(selectiontext2) 

    else:
        sel1 = bondvec[0]
        sel2 = bondvec[1]
        # if bond vectors are already specified with integer arrays
        try:
            bondvec_ndx = np.zeros([len(sel1), 2], dtype=np.int)
            bondvec_ndx[:,0] = np.array(sel1)
            bondvec_ndx[:,1] = np.array(sel2)
            return bondvec_ndx
        # if bond vectors are specified with selection texts
        except:
            atom1_ndx        = trj.topology.select(sel1)
            atom2_ndx        = trj.topology.select(sel2)  

    bondvec_ndx      = np.zeros([atom1_ndx.shape[0], 2], dtype=np.int)
    bondvec_ndx[:,0] = atom1_ndx
    bondvec_ndx[:,1] = atom2_ndx

    return bondvec_ndx


# ============================================================================ #


def _bond_vec(trj, bondvec_ndx):
    """
    Compute bond vectors from trajectory.

    Parameters
    ----------
    trj : mdtraj trajectory

    bondvec_ndx : (nvec, 2) array
        Bond vector indices

    Returns
    -------
    bondvec : (nframes, nvec, 3) array
        Array of normalized bond vectors
    info : dict
        Dictionary with information on the bondvectors.
        Keys are:
        'dt'            time step in pico seconds
        'bondlength'    one bondlength per bond vector
        'resnames'      two sequences of resnames, one for each atom in bondvector
        'resids'                         resids
        'resindex'                       resindices
        'atomnames'                      atomnames
        'atomindex'                      atomindices
        'element'                        elements
        'chain'                          chains ids
    """

    # extract bond vector subtrajectories
    atom1_xyz = trj.xyz[:,bondvec_ndx[:,0],:]
    atom2_xyz = trj.xyz[:,bondvec_ndx[:,1],:]

    # compute bond vectors
    bondvec = atom1_xyz - atom2_xyz

    # normalize bond vectors
    bondlength = (((bondvec**2).sum(2))**0.5)
    bondvec /= bondlength.reshape([bondvec.shape[0], bondvec.shape[1], 1])

    # estimate S2 with ensemble average formula from:
    # Trbovic et al. Proteins (2008). doi:10.1002/prot.21750
    S2 = np.zeros_like(bondvec[0,:,0])
    for i in range(3):
        for j in range(3):
            S2 += (bondvec[:,:,i] * bondvec[:,:,j]).mean(0)**2
    S2 = list(0.5 * (3*S2 - 1))

    info = {}
    info['dt'        ] = trj.timestep
    info['S2'        ] = S2
    info['bondlength'] = bondlength
    info['resnames'  ] = [[], []]
    info['resid'     ] = [[], []]
    info['resindex'  ] = [[], []]
    info['atomnames' ] = [[], []]
    info['atomindex' ] = [[], []]
    info['element'   ] = [[], []]
    info['chain'     ] = [[], []]

    # get info on atoms and residues
    for atom1_ndx, atom2_ndx in bondvec_ndx:
        atom1 = trj.topology.atom(atom1_ndx)
        atom2 = trj.topology.atom(atom2_ndx)
        info['resnames'  ][0].append(atom1.residue.name)
        info['resnames'  ][1].append(atom2.residue.name)
        info['resid'     ][0].append(atom1.residue.resSeq)
        info['resid'     ][1].append(atom2.residue.resSeq)
        info['resindex'  ][0].append(atom1.residue.index)
        info['resindex'  ][1].append(atom2.residue.index)
        info['atomnames' ][0].append(atom1.name)
        info['atomnames' ][1].append(atom2.name)
        info['atomindex' ][0].append(atom1.index) 
        info['atomindex' ][1].append(atom2.index) 
        info['element'   ][0].append(atom1.element) 
        info['element'   ][1].append(atom2.element)  
        info['chain'     ][0].append(atom1.residue.chain.index) 
        info['chain'     ][1].append(atom2.residue.chain.index)   

    return bondvec, info
 

# ============================================================================ #

 
def _angular_correlation_function(bondvec, verbose=True):
    """
    Compute 2nd order Legandre polynomial correlation functions of bond vectors.

    Parameters
    ----------
    bondvec : (nframes, nvectors, ndim) array
        normalized bond vectors
    verbose : boolean, optional
        Print progress messages or not

    Returns
    -------
    corr : (nvectors, nframes/2) array
        correlation functions
    corr_std : (nvectors, nframes/2) array
        standard deviations of correlation values
    corr_stdmean : (nvectors, nframes/2) array
        standard error of the mean
    """

    # compute bond vector correlation functions
    nframes, nvectors, ndim  = bondvec.shape
    corrnframes   = int(round(0.5*nframes))
    corr          = np.zeros([nvectors, corrnframes]) # correlation functions
    corr_std      = np.zeros([nvectors, corrnframes]) # standard deviation
    corr_stdmean  = np.zeros([nvectors, corrnframes]) # standard error of the mean
 
    for vector in range(nvectors):

        if verbose:
            print("\rProgress: {:3.0f}%".format(100.0*(vector+1)/nvectors), end="")
            sys.stdout.flush()

        # second order legendre polynomial of NH vector dotproducts P2[i,j] <=> P2[t,t+tau]
        P2 = np.polynomial.legendre.legval(np.dot(bondvec[:,vector,:], bondvec[:,vector,:].T), [0,0,1])

        # compute the correlation function for each lag time
        for frame in range(corrnframes):
            d = np.diagonal(P2, frame)
            corr        [vector, frame] = d.mean()
            corr_std    [vector, frame] = d.std()
            corr_stdmean[vector, frame] = corr_std[vector, frame] / d.shape[0]**0.5 

    if verbose:
        print()

    return corr, corr_std, corr_stdmean

 
# ============================================================================ #


def save_corr(savefilename, corr, corrstd, corrstdmean, bondvecinfo, topfilename, trjfilename, frames):
    """
    Save correlation functions to disk including all additional information.

    Parameters
    ----------
    savefilename : string
        Filename of the resulting zipfile containing all the data

    corr : (nvec, nframes) array
        Angular correlation functions for each bond vector

    corrstd : (nvec, nframes) array
        Angular correlation function standard deviations

    corrstdmean : (nvec, nframes) array
        Angular correlation function standard error of the mean

    bondvecinfo : dict
        Dictionary with information on bondlength, resids and resnames

    trjfilename : string
        Filename of the trajectory the original MD data was taken from.
        Stored as is, so better provide a complete path.

    topfilename : string
        Filename of the topology the original MD data was taken from.
        Stored as is, so better provide a complete path.

    frames : (startframe, endframe) tuple
        Range of frames used for correlation function computation
    """

    tmpdirnamelen     = 8
    tmpdirname        = "tmpdir_" + ''.join(np.random.choice([i for i in string.ascii_letters + string.digits], size=tmpdirnamelen))
    os.mkdir(tmpdirname)
    savefilepath      = os.path.dirname(savefilename)
    npzfilename       = "corr.npz"
    infofilename      = "info.dat"

    if not os.path.isdir(savefilepath) and not savefilepath == '':
        os.mkdir(savefilepath)

    info = {}
    info['npzfilename'] = npzfilename
    info['bondvecinfo'] = bondvecinfo
    info['trjfilename'] = trjfilename
    info['topfilename'] = topfilename
    info['frames']      = frames
 
    # save arrays
    with open(tmpdirname + '/' + npzfilename, 'wb') as outfile:
        np.savez_compressed(outfile, corr=corr, corrstd=corrstd, corrstdmean=corrstdmean)

    # save info
    with open(tmpdirname + '/' + infofilename, 'wb') as outfile:
        outfile.write(bz2.compress(pickle.dumps(info)))

    # pack both into zipfile
    with zipfile.ZipFile(savefilename, 'w') as outfile:
        outfile.write(tmpdirname + '/' + infofilename, arcname=infofilename)
        outfile.write(tmpdirname + '/' +  npzfilename, arcname=npzfilename )

    # remove npz and pickle files
    os.remove(tmpdirname + '/' + infofilename)
    os.remove(tmpdirname + '/' +  npzfilename)
    os.rmdir (tmpdirname)


# ============================================================================ #


def load_corr(filename):
    """
    Load correlation functions including all additional information.

    Parameters
    ----------
    filename : string
        Filename of the zipfile containing all the data

    Returns
    -------
    corr : (nvec, nframes) array
        Angular correlation functions for each bond vector

    corrstd : (nvec, nframes) array
        Angular correlation function standard deviations

    corrstdmean : (nvec, nframes) array
        Angular correlation function standard error of the mean

    info
        Dictionary with information on the loaded correlation functions
    """

    # create and change to working directory
    cwd           = os.getcwd()
    absfilename   = os.path.abspath(filename)
    tmpdirnamelen = 8
    tmpdirname    = "/tmp/" + ''.join(np.random.choice([i for i in string.ascii_letters + string.digits], size=tmpdirnamelen))
    os.mkdir(tmpdirname)
    os.chdir(tmpdirname)

    # extract files
    with zipfile.ZipFile(absfilename, 'r') as infile:
        zipfilenames = infile.namelist()
        for zipfilename in zipfilenames:
            infile.extract(zipfilename)

    for zipfilename in zipfilenames:
        with open(zipfilename, 'rb') as infile:
            # load info dictionary
            if zipfilename.find('.dat') >= 0:
                info = pickle.loads(bz2.decompress(infile.read()))

            # load corr arrays
            if zipfilename.find('.npz') >= 0:
                arrays      = np.load(infile)
                corr        = arrays['corr']
                corrstd     = arrays['corrstd']
                corrstdmean = arrays['corrstdmean']

        # remove extracted files
        os.remove(zipfilename)

    # change back to original cwd
    os.chdir(cwd)
    os.rmdir(tmpdirname)

    return corr, corrstd, corrstdmean, info


# ============================================================================ #
 

def bondvec_corr(trj, bondvec=None, fitgroup=None, parallel=True, saveinfo=None, verbose=True):
    """
    Compute bond vector correlation functions from mdtraj trajectory.

    Parameters
    ----------
    trj : mdtraj trajectory
        Trajectory to compute the order parameters of.

    bondvec : 2 element list/tuple, optional
        Both elements specify atoms to compute the bond vector from, either as mdtraj selection
        strings or as arrays with integer atom indices. Defaults to NH protein backbone bond vectors.

    fitgroup : selection text or numpay array, optional
        Atoms for RMSD fit, either specified by a selection string understandable by mdtraj or
        by an array containing integer atom indices.
        If not specified no fitting is done.

    parallel : boolean, optional
        Do the fitting in parallel on as many processors as available on the machine.

    saveinfo : dict, optional
        Information needed to save the correlation functions.
        Required keys are:
            'topfilename'   the topology   file the data comes from
            'trjfilename'   the trajectory file the data comes from
            'frames'        tuple with initial and final frame
            'zipfilename'   name of the zipfile to store correlation functions in

    verbose : boolean, optional
        Toggle verbosity level.


    Returns
    -------
    corr : (nvec, t) array
        Bond vector correlation functions.

    bondvecinfo :
        Dictionary with information on the bond vectors.
    """

    # fit trajectory
    if fitgroup is not None:
        _fit_trj(trj, fitgroup=fitgroup, parallel=parallel)

    bondvec_ndx                = _vec_atoms(trj, bondvec)
    bondvectors, bondvecinfo   = _bond_vec(trj, bondvec_ndx)
    corr, corrstd, corrstdmean = _angular_correlation_function(bondvectors, verbose=verbose)

    # store additional bond vector information
    if type(fitgroup) != type(None):
        bondvecinfo['fit' ] = True
    else:
        bondvecinfo['fit' ] = False
    bondvecinfo['fitgroup'] = fitgroup
    bondvecinfo['bondvec' ] = bondvec

    if type(saveinfo) == type(dict()):
        save_corr(saveinfo['zipfilename'], corr, corrstd, corrstdmean, bondvecinfo, saveinfo['topfilename'], saveinfo['trjfilename'], saveinfo['frames'])

    return corr, bondvecinfo

    
# ============================================================================ #

def bondvec_corr_batch_mpi(topfilename, trjfilenames, savepath, subtrjlength=None, bondvec=None, fitgroup=None):
    """
    Compute bond vector correlation functions from mdtraj trajectory.
    MPI parallelized batch version that stores the results on disk.

    Only present for compatibility with old scripts.
    Use bondvec_corr_batch() instead.
    """
 
    bondvec_corr_batch(topfilename, trjfilenames, savepath, subtrjlength=None, bondvec=None, fitgroup=None)

# ============================================================================ #

def bondvec_corr_batch(topfilename, trjfilenames, savepath, subtrjlength=None, bondvec=None, fitgroup=None):
    """
    Compute bond vector correlation functions from mdtraj trajectory.
    MPI parallelized batch version that stores the results on disk.

    Parameters
    ----------
    topfilename : string
        topology filename, e.g. PDB

    trjfilenames : sequence of strings
        trajectory filenames, e.g. XTC

    savepath : string
        path to store the results at

    subtrjlength: float, optional
        Subtrajetory length in ps.
        If none, whole trajectories are used.

    bondvec : 2 element list/tuple, optional
        Both elements specify atoms to compute the bond vector from, either as mdtraj selection
        strings or as arrays with integer atom indices. Defaults to NH protein backbone bond vectors.

    fitgroup : selection text or numpay array, optional
        Atoms for RMSD fit, either specified by a selection string understandable by mdtraj or
        by an array containing integer atom indices.
        If not specified no fitting is done.
    """

    # set up some timer and counter
    tc = {}
    tc['runtimer']  = time.time()
    tc['loadtimer'] = 0.0
    tc['corrtimer'] = 0.0
    tc['nsubtrjs' ] = 0
    tc['nframes'  ] = 0
 

    # set up MPI if present
    root   = 0
    if gotMPI:
        comm   = MPI.COMM_WORLD
        nranks = comm.Get_size()
        myrank = comm.Get_rank()
    else:
        comm   = None
        nranks = 1
        myrank = 0

    # to fit or not to fit?
    fit = "nofit"
    if type(fitgroup) != type(None):
        fit = "fit"

    # distribute tasks over ranks
    if myrank == root:
        ntasks        = len(trjfilenames)
        ntasksperrank = np.array(nranks * [ntasks / nranks], dtype=np.int)
        ntasksperrank[:ntasks%nranks] += 1

        for rank in range(1, nranks):
            taskstaken = sum(ntasksperrank[:rank])
            task = {}
            task['topfilename']  = topfilename
            task['trjindices']   = range(rank, len(trjfilenames), nranks)
            task['trjfilenames'] = [trjfilenames[i] for i in task['trjindices']]
            task['savepath']     = savepath

            comm.send(task, dest=rank, tag=rank)

    if myrank != root:
        task = comm.recv(source=root, tag=myrank)
    else:
        rank = myrank
        taskstaken = sum(ntasksperrank[:rank])
        task = {}
        task['topfilename']  = topfilename
        task['trjindices']   = range(rank, len(trjfilenames), nranks)
        task['trjfilenames'] = [trjfilenames[i] for i in task['trjindices']]
        task['savepath']     = savepath 


    # do the assigned piece of work
    for nf, trjfilename in enumerate(task['trjfilenames']):
        # determinde dt and chunksize
        trjs      = md.iterload(trjfilename, top=task['topfilename'], chunk=2)
        trj       = next(trjs)
        dt        = trj.timestep
        chunksize = int(subtrjlength / dt)
        trjindex  = task['trjindices'][nf]

        if myrank == root:
            print("\r", 100*" ", end="")
            print("\rProgress overall: {:3.0f}%; Current trajectory: ".format(100.0*nf/len(task['trjfilenames']), nf), end="")
            sys.stdout.flush()
 
        loadstarttime = time.time()
        for ntrj, trj in enumerate(md.iterload(trjfilename, top=task['topfilename'], chunk=chunksize)):

            tc['loadtimer'] += time.time() - loadstarttime

            if trj.n_frames == chunksize:
                tc['nsubtrjs' ] += 1
                tc['nframes']   += trj.n_frames

                saveinfo = {}
                saveinfo['topfilename'] = task['topfilename']
                saveinfo['trjfilename'] = trjfilename
                saveinfo['frames'     ] = (ntrj * chunksize, (ntrj+1)*chunksize-1)
                saveinfo['zipfilename'] = savepath + '/' + ''.join(os.path.basename(trjfilename).split('.')[:-1]) + '_{}.zip'.format(fit)
                saveinfo['zipfilename'] = "{}/{}_no{}_frames{}to{}_{}.zip".format(savepath,
                                                                                  ''.join(os.path.basename(trjfilename).split('.')[:-1]),
                                                                                  trjindex,
                                                                                  saveinfo['frames'][0],
                                                                                  saveinfo['frames'][1],
                                                                                  fit)
            
                corrstarttime = time.time()
                bondvec_corr(trj, bondvec=bondvec, fitgroup=fitgroup, parallel=False, saveinfo=saveinfo, verbose=False)
                tc['corrtimer'] += time.time() - corrstarttime
            loadstarttime = time.time()

            if myrank == root:
                print(".", end="")
                sys.stdout.flush()

    if myrank == root:
        print("")


    # report timers and counters to root
    if myrank == root:
        for rank in range(1, nranks):
            tci = comm.recv(source=rank, tag=rank)
            tc['loadtimer'] += tci['loadtimer'] 
            tc['corrtimer'] += tci['corrtimer'] 
            tc['nsubtrjs' ] += tci['nsubtrjs' ] 
            tc['nframes'  ] += tci['nframes'  ] 
    else:
        comm.send(tc, dest=root, tag=myrank)

    if gotMPI:
        comm.barrier()

    tc['runtimer'] = time.time() - tc['runtimer']
    if myrank == root:
        print("Loaded {} subtrajectories with a total of {} frames".format(tc['nsubtrjs'], tc['nframes']))
        print("Total runtime:                        {:8.0f} sec.".format(tc['runtimer' ]))
        

