# -*- coding: utf-8 -*-
# @Author: Helios
# @Date:   2017-07-13 14:20:04
# @Last Modified by:   Helios
# @Last Modified time: 2017-10-10 22:06:20


# this entire script is just gross but almost all of the functionality was
# meant to be one off use, hence...this.

import os
import sys
import h5py
import matplotlib
import numpy as np
import matlab.engine
from time import gmtime
import ibmtomography as itm
from itertools import product
import matplotlib.pyplot as plt


# path to experiment data file
if os.path.isdir('C:\\Users\\Helios'):
    archivepath = 'C:\\Users\\Helios\\Documents\\Projects\\Uni\\2017\\markovpy\\archive.hdf5'
else:
    archivepath = 'C:\\Users\\Joshua\\Documents\\Projects\\Uni\\2017\\Research\\markovpy\\archive.hdf5'


#--------------------------------------------------------------------------
# COMPLETE
#--------------------------------------------------------------------------


class QCircExp(object):
    """
    Defines a custom data structure to store specific quantum circuits and associated 
    reconstruction information
    """

    def __init__(self, archive, circuit):
        # store circuit name
        self.name = circuit
        # retrieve and store circuit attributes
        self.device = archive[circuit].attrs['device'].decode('ascii')
        self.tomography = archive[circuit].attrs['tomography'].decode('ascii')
        self.qubits = archive[circuit].attrs['qubits'][()]
        self.complete = archive[circuit].attrs['complete'][()]
        self.total = archive[circuit].attrs['total'][()]
        self.rawqasm = archive[circuit]['raw_qasm'][()]
        self.shots = archive[circuit].attrs['shots'][()]

        # check for circuit reconstruction data structures
        if self.complete == self.total:
            # store key information from either
            if self.tomography == 'state':
                self.finalstate = archive[circuit]['tomography_ML'][()]
            elif self.tomography == 'process':
                # extract Chi matrix
                self.chi = np.asarray(archive[circuit]['Data_Group'][
                                      'Process_matrix'][()])
                # Kraus operators
                self.kraus = np.asarray(
                    archive[circuit]['Data_Group']['Kraus_set'][()])
                # operator basis used in reconstruction
                self.opbasis = np.asarray(archive[circuit]['Data_Group'][
                                          'Operator_basis_set'][()])
                # preperation basis used in reconstruction
                self.prepbasis = np.asarray(archive[circuit]['Data_Group'][
                                            'Preperation_basis_set'][()])
                # Choi matrix
                self.choi = np.asarray(
                    archive[circuit]['Data_Group']['Choi_matrix'][()])
                # extract tomography dataset
                self.tomographyset = []
                for state in archive[circuit]['tomography_ML']:
                    self.tomographyset.append(
                        archive[circuit]['tomography_ML'][state][()])


# imports target data from archive and converts to analysis class
def circuitimport(circuit):
    with h5py.File(archivepath, 'a') as archive:
        # create shortcut variable and initialise class
        return QCircExp(archive, circuit)

# compute trace distance of two matrices using good ol trace norm


def tracedistcrude(m1, m2):
    dist = 0.5*np.trace(np.sqrt(np.asmatrix(m1-m2).H*(m1-m2)))
    return np.real(dist)


# compares actual and simulator results process tomography using some base name
def comparechiplot(basename, circuits=10, method='diamond'):
    # iterate over different versions and add circuit information to primary
    # data container
    circcollect = []
    for i in range(1, circuits+1):
        # simulator path
        simname = basename + '{}'.format(i) + '_simulator'
        # equivalent hardware path
        hardname = basename + '{}'.format(i) + '_ibmqx2'
        # construct circuit classes
        simexp = circuitimport(simname)
        hardexp = circuitimport(hardname)
        # store as immutable type
        circcollect.append(tuple([simexp, hardexp]))
    circcollect = np.asarray(circcollect)

    # compute difference vector using input method
    channeldist = np.asarray([])
    if method == 'trace':
        for pair in circcollect:
            channeldist = np.append(
                channeldist, tracedistcrude(pair[0].chi, pair[1].chi))
    elif method == 'diamond':
        # import and initialise a matlab engine instance
        import matlab.engine
        mateng = matlab.engine.start_matlab()
        # compute diamond norm for 'equivalent' quantum channels
        for pair in circcollect:
            diffchoi = pair[0].choi - pair[1].choi
            diffchoi = matlab.double(diffchoi.tolist(), is_complex=True)
            channeldist = np.append(channeldist, mateng.dnorm2(diffchoi))
    else:
        print('------- UNKNOWN COMPARISON METHOD {}: ABORTING -------'.format(method))
        exit()

    # compute crude statistical bound lines
    av = np.mean(channeldist)
    statdev = 1/np.sqrt(circcollect[0][1].shots)
    # compute crudely fitted bound lines
    gatenum = list(range(1, len(channeldist)+1))
    bnds = np.poly1d(np.polyfit(gatenum, channeldist, 2))
    print(statdev)
    upbnd = bnds(gatenum) + statdev
    lowbnd = bnds(gatenum) - statdev

    plt.plot(channeldist, label='Distance')
    plt.plot(upbnd, 'r--', label='Bounds')
    plt.plot(lowbnd, 'r--')
    plt.ylim([0, 0.1])
    plt.title(
        'Channel distance as a function of gate number for circuit: ' + basename)
    plt.grid(True)
    plt.xlabel('Number of applied random unitaries')
    plt.ylabel('Channel distance using ' + method + ' norm')
    plt.legend()
    plt.show()


def tracenormcompare(basename, circuits=10, method='sim'):
    # iterate over different versions and add circuit information to primary
    # data container
    circcollect = []
    for i in range(1, circuits+1):
        # simulator path
        simname = basename + '{}'.format(i) + '_simulator'
        # equivalent hardware path
        hardname = basename + '{}'.format(i) + '_ibmqx2'
        # construct circuit classes
        simexp = circuitimport(simname)

        hardexp = circuitimport(hardname)
        # store as immutable type
        circcollect.append(tuple([simexp, hardexp]))
    circcollect = np.asarray(circcollect)

    # compare trace norm of density matricies for two states
    tracedist = np.asarray([])
    if method == 'sim':
        for pair in circcollect:
            #print(tracenorm(pair[0].tomographyset[0] - pair[0].tomographyset[1]))
            tracedist = np.append(tracedist, tracenorm(
                pair[0].tomographyset[0] - pair[0].tomographyset[1]))
    else:
        for pair in circcollect:
            tracedist = np.append(tracedist, tracenorm(
                pair[1].tomographyset[0] - pair[1].tomographyset[1]))

    gatenum = list(range(1, len(tracedist)+1))
    plt.plot(gatenum, tracedist)
    plt.title(
        'Trace norm as a function of hadamard gates applied for circuit: ' + basename)
    plt.grid(True)
    plt.ylim([1.8, 2.1])
    plt.xlabel('Number of applied gates')
    plt.ylabel('distance using tracenorm')
    plt.legend('Actual', 'Ideal')
    plt.show()


def tracenormcomparetmp(circuits=10, method='sim'):
    # iterate over different versions and add circuit information to primary
    # data container
    circcollect = []
    for i in range(1, circuits+1):
        # hardware path
        hardname1 = '1QPT_Hadamard1_H' + '{}'.format(i) + '_ibmqx2'
        hardname2 = '1QPT_Hadamard2_H' + '{}'.format(i) + '_ibmqx2'
        hardname3 = '1QPT_Hadamard3_H' + '{}'.format(i) + '_ibmqx2'
        hardname4 = '1QPT_Hadamard4_H' + '{}'.format(i) + '_ibmqx2'

        hardexp1 = circuitimport(hardname1)
        hardexp2 = circuitimport(hardname2)
        hardexp3 = circuitimport(hardname3)
        hardexp4 = circuitimport(hardname4)

        # store as immutable type
        circcollect.append(tuple([hardexp1, hardexp2, hardexp3, hardexp4]))

    circcollect = np.asarray(circcollect)

    # compare trace norm of density matricies for two states
    tracedist = np.asarray([])

    for group in circcollect:
        average0 = group[0].tomographyset[0] + group[1].tomographyset[0] + \
            group[2].tomographyset[0] + group[3].tomographyset[0]
        average0 = average0/4.0
        average1 = group[0].tomographyset[1] + group[1].tomographyset[1] + \
            group[2].tomographyset[1] + group[3].tomographyset[1]
        average1 = average1/4.0

        tracedist = np.append(tracedist, tracenorm(average0 - average1))

    gatenum = list(range(1, len(tracedist)+1))
    plt.plot(gatenum, tracedist, label='Actual', linewidth=3.0)
    plt.title(
        'Trace norm as a function of hadamard gates applied for circuit: 1QPT_Hadamard')
    plt.grid(True)
    plt.ylim([1.8, 2.05])
    plt.plot(gatenum, [2.0]*len(tracedist),
             label='Ideal', linewidth=3.0, ls='dashed')
    plt.xlabel('Number of applied gates')
    plt.ylabel('$|| \Lambda(|0><0|) - \Lambda(|1><1|)||_1$')
    matplotlib.rcParams.update({'font.size': 22})
    plt.legend()
    plt.show()


def rhoplot(rho, axislabels=None, save=False):
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    from mpl_toolkits.mplot3d import Axes3D

    # extract real and imaginary components of density matrix
    realrho = np.real(rho)
    imagrho = np.imag(rho)

    # instantiate new figure
    fig = plt.gcf()
    fig.canvas.set_window_title('Density Plot')
    #rax = Axes3D(fig)
    rax = fig.add_subplot(121, projection='3d')
    iax = fig.add_subplot(122, projection='3d')

    # set titles
    rax.title.set_text('Real$(\\rho)$')
    iax.title.set_text('Imag$(\\rho)$')
    # apply custom labelling
    if axislabels is not None:
        rax.set_xticklabels(axislabels)
        rax.set_yticklabels(axislabels)
        iax.set_xticklabels(axislabels)
        iax.set_yticklabels(axislabels)

    # dimension of space
    dim = np.shape(realrho)[0]
    # create indexing vectors
    x, y = np.meshgrid(range(0, dim), range(0, dim), indexing='ij')
    x = x.flatten('F')
    y = y.flatten('F')
    z = np.zeros_like(x)

    # create bar widths
    dx = 0.5*np.ones_like(z)
    dy = dx.copy()
    dzr = realrho.flatten()
    dzi = imagrho.flatten()

    # compute colour matrix for real matrix and set axes bounds
    norm = colors.Normalize(dzr.min(), dzr.max())
    rcolours = cm.Reds(norm(dzr))
    rax.set_zlim3d([0, np.max(dzr)])
    iax.set_zlim3d([0, np.max(dzr)])

    inorm = colors.Normalize(dzi.min(), dzi.max())
    icolours = cm.jet(inorm(dzi))

    # plot image
    rax.bar3d(x, y, z, dx, dy, dzr, color=rcolours)
    iax.bar3d(x, y, z, dx, dy, dzi, color=icolours)
    #plt.ticklabel_format(style='sci', axis='z', scilimits=(0, 0))
    plt.show()


def condplot(errors, axislabels=None, save=False):
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    from mpl_toolkits.mplot3d import Axes3D

    # extract real and imaginary components of density matrix

    # instantiate new figure
    fig = plt.gcf()
    fig.canvas.set_window_title('Density Plot')
    #rax = Axes3D(fig)
    rax = fig.add_subplot(111, projection='3d')

    # set titles
    rax.title.set_text('Conditional Error Rates')
    # apply custom labelling
    if axislabels is not None:
        rax.set_xticklabels(axislabels)
        rax.set_yticklabels(axislabels)

    # dimension of space
    dim = np.shape(errors)[0]
    # create indexing vectors
    x, y = np.meshgrid(range(0, dim), range(0, dim), indexing='ij')
    x = x.flatten('F')
    y = y.flatten('F')
    z = np.zeros_like(x)

    # create bar widths
    dx = 0.5*np.ones_like(z)
    dy = dx.copy()
    dzr = errors.flatten()

    # compute colour matrix for real matrix and set axes bounds
    norm = colors.Normalize(dzr.min(), dzr.max())
    rcolours = cm.Reds(norm(dzr))
    rax.set_zlim3d([0, np.max(dzr)])
    rax.set_xlabel('$V$', fontsize=20)
    rax.set_ylabel('$U$', fontsize=20)

    # plot image
    rax.bar3d(x, y, z, dx, dy, dzr, color=rcolours)
    #plt.ticklabel_format(style='sci', axis='z', scilimits=(0, 0))
    plt.show()


def densityplot(archive, circuit, method="ML", save=False):
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    from mpl_toolkits.mplot3d import Axes3D

    rho = archive[circuit]['tomography_' + method]

    # extract real and imaginary components of density matrix
    realrho = np.real(rho)
    imagrho = np.imag(rho)

    # instantiate new figure
    fig = plt.gcf()
    fig.canvas.set_window_title('Density Plot for circuit {}'.format(circuit))
    #rax = Axes3D(fig)
    rax = fig.add_subplot(121, projection='3d')
    iax = fig.add_subplot(122, projection='3d')

    # set titles
    rax.title.set_text('Real$(\\rho)$')
    iax.title.set_text('Imag$(\\rho)$')

    # dimension of space
    dim = np.shape(realrho)[0]
    # create indexing vectors
    x, y = np.meshgrid(range(0, dim), range(0, dim), indexing='ij')
    x = x.flatten('F')
    y = y.flatten('F')
    z = np.zeros_like(x)

    # create bar widths
    dx = 0.5*np.ones_like(z)
    dy = dx.copy()
    dzr = realrho.flatten()
    dzi = imagrho.flatten()

    # compute colour matrix for real matrix and set axes bounds
    norm = colors.Normalize(dzr.min(), dzr.max())
    rcolours = cm.Reds(norm(dzr))
    rax.set_zlim3d([0, 0.25])
    iax.set_zlim3d([0, 0.25])

    inorm = colors.Normalize(dzi.min(), dzi.max())
    icolours = cm.jet(inorm(dzi))

    # plot image
    rax.bar3d(x, y, z, dx, dy, dzr, color=rcolours)
    iax.bar3d(x, y, z, dx, dy, dzi, color=icolours)
    plt.ticklabel_format(style='sci', axis='z', scilimits=(0, 0))
    plt.show()

# IBM arbitrary unitary gate


def ibmu3(theta, phi, lmbda):
    return np.asarray([[np.cos(theta/2), -np.exp(1j*lmbda)*np.sin(theta/2)],
                       [np.exp(1j*phi)*np.sin(theta/2), np.exp(1j*lmbda + 1j*phi)*np.cos(theta/2)]])


def chiplot(archive, circuit):
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    from mpl_toolkits.mplot3d import Axes3D

    chi = archive[circuit]['Data_Group']['Process_matrix'][()]

    # extract real and imaginary components of density matrix
    realchi = np.real(chi)
    imagchi = np.imag(chi)

    # instantiate new figure
    fig = plt.gcf()
    fig.canvas.set_window_title(
        'Process matrix for circuit {}'.format(circuit))
    #rax = Axes3D(fig)
    rax = fig.add_subplot(121, projection='3d')
    iax = fig.add_subplot(122, projection='3d')

    # set titles
    rax.title.set_text('Real$(\\chi)$')
    iax.title.set_text('Imag$(\\chi)$')

    # dimension of space
    dim = np.shape(realchi)[0]
    # create indexing vectors
    x, y = np.meshgrid(range(0, dim), range(0, dim), indexing='ij')
    x = x.flatten('F')
    y = y.flatten('F')
    z = np.zeros_like(x)

    # create bar widths
    dx = 0.5*np.ones_like(z)
    dy = dx.copy()
    dzr = np.abs(realchi.flatten())
    dzi = np.abs(imagchi.flatten())

    # compute colour matrix for real matrix and set axes bounds
    norm = colors.Normalize(dzr.min(), dzr.max())
    rcolours = cm.Reds(norm(dzr))
    rax.set_zlim3d([0, dzr.max()])
    iax.set_zlim3d([0, dzr.max()])

    inorm = colors.Normalize(dzi.min(), dzi.max())
    icolours = cm.jet(inorm(dzi))

    # plot image
    rax.bar3d(x, y, z, dx, dy, dzr, color=rcolours)
    iax.bar3d(x, y, z, dx, dy, dzi, color=icolours)
    plt.ticklabel_format(style='sci', axis='z', scilimits=(0, 0))
    plt.show()


def tracenorm(m):
    import numpy.linalg
    return np.sum(np.abs(numpy.linalg.eigh(m)[0]))


# compute A form of Choi state (also computes the choi state if given the
# A form)
def choi2A(choi):
    # retrieve dimension
    dim = len(choi)
    # reshape dimension
    sdim = int(np.sqrt(dim))
    loui = np.asarray(choi)
    return np.transpose(loui.reshape([sdim, sdim, sdim, sdim]), [0, 2, 1, 3]).reshape([dim, dim])


# compute the conditional diamond distance of the various gate combinations
def conditional(archive):
    # initialise matlab engine
    mateng = matlab.engine.start_matlab()
    # define basis set
    gateset = ['X', 'Y', 'Z', 'T', 'H', 'S', 'CX']
    #gateset = ['X','Z','CX']
    # construct all combinations of the above
    #combs = product(gateset, repeat=2)

    errormatrix = np.zeros([len(gateset)]*2)
    for num1, iter1 in enumerate(gateset):
        for num2, iter2 in enumerate(gateset):
            try:
                errormatrix[num1][num2] = condell(
                    archive, [iter1, iter2], mateng)
            # catch incomplete combinations
            except KeyError:
                errormatrix[num1][num2] = 0.0
    return errormatrix, gateset


# Computes the conditional distance of an operator
def condell(archive, item, mateng):
    # extract gate combination
    initgate = item[0]
    fingate = item[1]

    # generate circuit paths (this entire function is so shit)
    if initgate == 'CX':
        ibmcomb = '2QPT' + '_' + initgate + fingate + '_ibmqx4'
        ibminit = '2QPT' + '_' + initgate + '_ibmqx4'
        simfin = '1QPT' + '_' + fingate + '_simulator'
        twoflag = 1
    elif fingate == 'CX':
        ibmcomb = '2QPT' + '_' + initgate + fingate + '_ibmqx4'
        simfin = '2QPT' + '_' + fingate + '_simulator'
        ibminit = '1QPT' + '_' + initgate + '_ibmqx4'
        twoflag = 2
    else:
        ibmcomb = '1QPT' + '_' + initgate + fingate + '_ibmqx4'
        simfin = '1QPT' + '_' + fingate + '_simulator'
        ibminit = '1QPT' + '_' + initgate + '_ibmqx4'
        twoflag = 0
    # generate everything from the Kraus set because mapping A form to higher
    # dimensional A form is hard
    if twoflag == 1:
        simchoifin = itm.numpyarr2matlab(itm.kraus2choi(
            archive[simfin]['Data_Group']['Kraus_set'][()], targets=2))
        initAhard = itm.kraus2choi(archive[ibminit]['Data_Group'][
                                   'Kraus_set'][()], rep='loui')
    elif twoflag == 2:
        simchoifin = itm.numpyarr2matlab(itm.kraus2choi(
            archive[simfin]['Data_Group']['Kraus_set'][()]))
        initAhard = itm.kraus2choi(archive[ibminit]['Data_Group'][
                                   'Kraus_set'][()], rep='loui', targets=2)
    else:
        initAhard = itm.kraus2choi(archive[ibminit]['Data_Group'][
                                   'Kraus_set'][()], rep='loui')
        simchoifin = itm.numpyarr2matlab(itm.kraus2choi(
            archive[simfin]['Data_Group']['Kraus_set'][()]))

    comboAhard = itm.kraus2choi(archive[ibmcomb]['Data_Group'][
                                'Kraus_set'][()], rep='loui')

    # compute inverse of the A form of initial map
    B = np.asarray(mateng.matrixinverse(itm.numpyarr2matlab(initAhard)))
    # compute actual distance of gate
    actual = itm.numpyarr2matlab(choi2A(np.dot(comboAhard, B)))
    dist = float(mateng.dnorm(actual, simchoifin))
    return dist



# flattens a list of lists into a depth one list
def listflatten(lists):
    return [item for sublist in lists for item in sublist]


# computes distance of a given process tensor from that of the closest markovian process
# using the chosen distance metric: quantum relative entropy or tracenorm
# (only included for thoroughness)
def markovdist(ptensor, subsystems, metric='qre', gen=False):
    # compute dimensionality of system (corresponds to number of timesteps)
    mptensor = 1.0
    # compute perms for later as partial trace deletes elements in subsys -
    # should I change that?
    perms = listflatten(subsystems)
    for item in subsystems:
        mptensor = np.kron(mptensor, itm.partialtrace(ptensor, item))
    
    # permute B form to match reconstructed markovian tensor
    ptensor = subsyspermute(ptensor, perms, [2, 2, 2, 2])
    #matrix2dat(ptensor, 'permutedsim.dat')
    # compute distance using chosen metric
    if metric == 'qre':
        # compute quantum relative entrop
        from scipy.linalg import logm
        dist = np.trace(np.dot(ptensor, (logm(ptensor)/np.log(2) - logm(mptensor)/np.log(2))))
        if gen:
            return dist, mptensor
        else:
            return dist
    elif metric == 'trnm':
        # compute trace norm
        return tracenorm(ptensor, mptensor)
    else:
        print('Unknown metric')
        return np.inf

# permutes the subsystems of a density matrix with dimensions dims
# according to perm


def subsyspermute(rho, perm, dims):
    # get dimensions of system
    d = np.shape(rho)
    # get number of subsystems
    sys = len(dims)
    # perform permutation
    perm = [(sys - 1 - i) for i in perm[-1::-1]]
    perm = listflatten([perm, [sys + j for j in perm]])
    return np.transpose(rho.reshape(dims[-1::-1]*2), perm).reshape(d)


# outputs a numpy array to a dat file accessible by mma. 
def matrix2dat(array, filename):
    array = np.asarray(np.asmatrix(array).flatten())
    array.astype('complex128').tofile(filename)


# computes the kraus operators given B form of a map
def choi2kraus(choi):
    # import singular value decompoisiton library
    from scipy.linalg import svd
    # retrieve eigenvectors/eigenvalues
    U,s,Vh = svd(choi)
    # get system dimension
    dim = int(np.log2(len(choi)))
    kraus = []
    for col,lam in enumerate(s):
        operator = np.sqrt(lam)*np.reshape(Vh[col][:],(dim,dim))
        kraus.append(operator)
    return kraus


# expands a single qubit operator to be applied to a larger space
def operatorexpand(op, target, total):
    newop = 1.0
    for qubit in range(0,total):
        if target == qubit:
            newop = np.kron(newop, op)
        else:
            newop = np.kron(newop, itm.gatedict['id'])
    return newop


# extends a quantum channel to include identity operations on the larger system
def choipad(choi, targets, total):
    # get total system dimensions
    dim = 2**total
    # get kraus operators for choi state
    kraus = choi2kraus(choi)
    # initialise 
    newkraus = []
    for op in kraus:
        newop = 1.0
        appflag = False
        for i in range(0, total):
            # slot in applied operator
            if i in targets:
                if not appflag:
                    appflag = True
                    newop = np.kron(newop, op)
            else:
                newop = np.kron(newop, itm.gatedict['id'])
        newkraus.append(newop)

    # compute corresponding B form
    newchoi = itm.kraus2choi(newkraus)
    return newchoi




#--------------------------------------------------------------------------
# PURGATORY
#--------------------------------------------------------------------------

# determines the CPTP maps that act on the subsystems of a perfect swap gate such that the resultant operation
# approximates the actual swap. Returns missing aspect in B form
def purifyswap(archive):
    # extract Choi states of Swap operation
    swapactual = archive['2QPT_SWAP_ibmqx2']['Data_Group']['Choi_matrix'][()]
    swapsim = archive['2QPT_SWAP_simulator']['Data_Group']['Choi_matrix'][()]

    # convert to A form and format for matlab
    swapactual = itm.numpyarr2matlab(choi2A(swapactual))
    swapsim = itm.numpyarr2matlab(choi2A(swapsim))
    swap = itm.numpyarr2matlab(itm.kraus2choi(np.asarray([[[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]]]), rep='loui'))

    # compute perturbation such that perturb*swap = swapactual
    mateng = matlab.engine.start_matlab()
    perturb = itm.numpyarr2matlab(np.asarray(mateng.liouestimate(swap, swapactual)))
    

    Tkraus = archive['1QPT_T_ibmqx4']['Data_Group']['Kraus_set'][()]
    Tkraus = [operatorexpand(i, 0, 2) for i in Tkraus]
    TAform = itm.kraus2choi(Tkraus, rep='loui')

    # fix process in A form
    fix = np.asarray(mateng.unperturb(perturb, itm.numpyarr2matlab(TAform)))

    #rhoplot(fix*np.dot(np.asarray(perturb),TAform))

    #matrix2dat(fix, 'swapfix.dat')

def processmap(archive, processtensor):
    # extract choi state of swap gate
    swapactual = itm.numpyarr2matlab(archive['2QPT_SWAP_ibmqx4']['Data_Group']['Choi_matrix'][()]) 
    # extract process tensor
    ptensorbad = itm.numpyarr2matlab(archive[processtensor]['tomography_ML'][()])
    # pad out choi state to act on 4 qubit state and convert to matlab type
    swapactualpad = itm.numpyarr2matlab(choipad(swapactual, [1,2], 4))
    #define perfect swap
    swap = itm.kraus2choi(np.asarray([[[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]]]), rep='choi')
    # initialise matlab engine
    mateng = matlab.engine.start_matlab()
    # compute map
    pmap = np.asarray(mateng.processmap(swapactual, ptensorbad))
    # apply map to perfect swap
    ptensor = mapapply(pmap, swap)
    # plot best guess at correct one
    #rhoplot(ptensorbad)
    print(np.trace(np.asarray(ptensor)))
    # compute QRE of new ptensor
    print(markovdist(np.asarray(ptensor)/np.trace(ptensor), [[1,2], [0,3]]))

    #rhoplot(np.asarray(ptensor)/np.trace(ptensor))
    rhoplot(itm.partialtrace(np.asarray(ptensor)/np.trace(ptensor), [1,2]))

def mapapply(pmap, control):
    # compute dimensionality
    dim = len(control)
    # vectorise control operation
    controlflat = np.reshape(control, [dim**2, 1])
    # apply map 
    outstate = np.dot(pmap, controlflat)
    # reshape and return
    return np.reshape(outstate, [dim, dim]).T

def choiexport(archive):
    #retrieve choi states of required dynamics
    phase       = archive['1QPT_S_ibmqx4']['Data_Group']['Choi_matrix'][()]
    phasesim    = archive['1QPT_S_simulator']['Data_Group']['Choi_matrix'][()]

    swap        = archive['2QPT_SWAP_ibmqx4']['Data_Group']['Choi_matrix'][()]
    swapsim     = archive['2QPT_SWAP_simulator']['Data_Group']['Choi_matrix'][()]

    tgate       = archive['1QPT_T_ibmqx4']['Data_Group']['Choi_matrix'][()]
    tgatesim    = archive['1QPT_T_simulator']['Data_Group']['Choi_matrix'][()]

    ptensor     = archive['ProcessTensorInternalSwap_ibmqx4']['tomography_ML'][()]
    ptensorsim  = archive['ProcessTensorInternalSwap_simulator']['tomography_ML'][()]

    ptensornorm     = archive['ProcessTensor_ibmqx4']['tomography_ML'][()]
    ptensornormsim  = archive['ProcessTensor_simulator']['tomography_ML'][()]

    cnot  = archive['2QPT_CX_ibmqx4']['Data_Group']['Choi_matrix'][()]
    cnotsim  = archive['2QPT_CX20_simulator']['Data_Group']['Choi_matrix'][()]

    had = archive['1QPT_H_ibmqx4']['Data_Group']['Choi_matrix'][()]
    hadsim = archive['1QPT_H_simulator']['Data_Group']['Choi_matrix'][()]

    pauliX = archive['1QPT_X_ibmqx4']['Data_Group']['Choi_matrix'][()]
    pauliXsim = archive['1QPT_X_simulator']['Data_Group']['Choi_matrix'][()]

    # store as dat files for mma access
    #matrix2dat(entangle,    'entangle.dat')
    #matrix2dat(entanglesim, 'entanglesim.dat')
    matrix2dat(phase,       'phase.dat')
    matrix2dat(phasesim,    'phasesim.dat')
    matrix2dat(swap,        'swap.dat')
    matrix2dat(swapsim,     'swapsim.dat')
    matrix2dat(tgate,       'tgate.dat')
    matrix2dat(tgatesim,    'tgatesim.dat')
    matrix2dat(ptensor,     'ptensor.dat')
    matrix2dat(ptensorsim,  'ptensorsim.dat')
    matrix2dat(ptensornorm,     'ptensornorm.dat')
    matrix2dat(ptensornormsim,  'ptensornormsim.dat')
    matrix2dat(cnot,        'cnot.dat')
    matrix2dat(cnotsim,     'cnotsim.dat')
    matrix2dat(had,         'had.dat')
    matrix2dat(hadsim,      'hadsim.dat')
    matrix2dat(pauliX,      'pauliX.dat')
    matrix2dat(pauliXsim,   'paulixsim.dat')







#--------------------------------------------------------------------------
# TO BE IMPLEMENTED
#--------------------------------------------------------------------------


# applies the map described by the process tensor
def choitensor(archive, circuit, rho):
    pass


# extracts the map associated with time step k in the Sudarshan B form of
# the process tensor
def processextract(archive, circuit, subsystem):
    pass


#--------------------------------------------------------------------------
# CURRENTLY WORKING ON
#--------------------------------------------------------------------------

def correlationplot(archive, circuit):
    ptensor = archive[circuit]['tomography_ML'][()]
    perms = [[0,2],[1,3]]
    mptensor = 1.0
    for item in perms:
        mptensor = np.kron(mptensor, itm.partialtrace(ptensor, item))

    corr = np.corrcoef(ptensor)
    mcorr = np.corrcoef(mptensor)
    rhoplot(corr - mcorr)
#--------------------------------------------------------------------------







if __name__ == '__main__':
    with h5py.File(archivepath, 'a') as archive:
        #correlationplot(archive, 'ProcessTensorInternalSwap_simulator')
        #choiexport(archive)
        #processmap(archive, 'ProcessTensor_ibmqx4')
        #purifyswap(archive)
        #mateng = matlab.engine.start_matlab()
        #array,t,w = itm.processtomography(archive, '1QPT_HH_simulator', mateng)
        #rhoplot(archive['Entangle_ibmqx4']['Data_Group']['Choi_matrix'][()] - archive['Entangle_simulator']['Data_Group']['Choi_matrix'][()])
        #matrix2dat(array, 'ASWAPsim.dat')
        #ptensor2 = choi2A(archive['ProcessTensorCorrelatedSwap5_simulator']['tomography_ML'][()])
        #ptensor4 = choi2A(archive['ProcessTensorCorrelatedSwap4_simulator']['tomography_ML'][()])
        #mateng = matlab.engine.start_matlab()
        #permute = np.asarray(mateng.permutecompute(itm.numpyarr2matlab(ptensor2),itm.numpyarr2matlab(ptensor4)))
        #iden = itm.partialtrace(ptensor, [1,2])
        #SWAP = archive['2QPT_SWAP_simulator']['Data_Group']['Process_matrix'][()]
        #SWAP = archive['2QPT_SWAP_simulator']['Data_Group']['Process_matrix'][()]
        #purifyswap(archive)
        #ptensor = archive['ProcessTensorXX_ibmqx4']['tomography_ML'][()]
        #rhoplot(ptensor)
        #print(markovdist(ptensor, [[1,2],[0,3]]))
        
        errors, gateset = conditional(archive)
        condplot(errors, gateset)
