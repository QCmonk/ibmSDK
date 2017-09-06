# -*- coding: utf-8 -*-
# @Author: Helios
# @Date:   2017-07-13 14:20:04
# @Last Modified by:   Helios
# @Last Modified time: 2017-09-06 20:10:19


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
        return QCircExp(archive[circuit], circuit)

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
            channeldist = np.append(channeldist, mateng.dnorm(diffchoi))
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


def rhoplot(rho, save=False):
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
    dzr = np.abs(realrho.flatten())
    dzi = np.abs(imagrho.flatten())

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
    plt.ticklabel_format(style='sci', axis='z', scilimits=(0, 0))
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
    dzr = np.abs(realrho.flatten())
    dzi = np.abs(imagrho.flatten())

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


#--------------------------------------------------------------------------
# TO BE IMPLEMENTED
#--------------------------------------------------------------------------

# compute the conditional diamond distance of the various gate combinations
def conditional(archive, base):
    # initialise matlab engine
    mateng = matlab.engine.start_matlab()
    # define basis set
    gateset = ['H', 'T', 'S', 'X', 'Y', 'Z']
    # construct all combinations of the above
    combs = product(gateset, repeat=2)
    for iteration in combs:
        pass


def condell(archive, base, item, mateng):
    # extract gate combination
    initgate = item[0]
    fingate = item[1]

    # retrieve required choi matrices
    combpathhard = base + '_' + initgate + fingate + '_ibmqx2'
    initpathsim = base + '_' + initgate + '_simulator'
    finpathhard = base + '_' + fingate + '_ibmqx2'

    combochoihard = archive[combpathhard]['Data_Group']['Choi_matrix'][()]
    initchoisim = archive[initpathsim]['Data_Group']['Choi_matrix'][()]
    finchoihard = itm.numpyarr2matlab(
        archive[finpathhard]['Data_Group']['Choi_matrix'][()])

    # compute inverse of final gate
    finchoihardinv = np.asarray(mateng.matrixinverse(finchoihard))

    diffchoi = itm.numpyarr2matlab(combochoihard*finchoihardinv - initchoisim)
    dist = mateng.dnorm(diffchoi)
    print(dist)

#--------------------------------------------------------------------------
# CURRENTLY WORKING ON
#--------------------------------------------------------------------------

# applies the map described by the process tensor


def choitensor(archive, circuit, rho):
    pass


# computes the partial trace of m \in \mathcal{H^n}, tracing out subsystems not in sys
# why must you make my life so difficult numpy?
def partialtrace(m, sys):
    # type enforcement
    m = np.asarray(m)
    # get tensor dimensions
    qnum = int(np.log2(len(m)))
    # compute dimensions of tensor
    tshape = (2,)*2*qnum
    # reshape to tensor
    mtensor = m.reshape((tshape))
    # compute dimensions to trace over
    index1, index2 = sys[0], sys[0] + qnum
    del sys[0]
    newdim = 2**(qnum-1)
    # compute reduced density matrix via recursion
    if len(sys) > 0:
        # trace out target subsystem (repeated reshaping is a bit rough but its not worth the fix time) 
        mtensor = np.trace(mtensor, axis1=index1, axis2=index2).reshape((newdim, newdim))
        # adjust subsequent target dimensions with shallow copy
        sys[:] = [i-1 for i in sys]
        # by the power of recursion
        mtensor = partialtrace(mtensor, sys)
    else:
        # bottom of the pile, compute and pass up the chain
        mtensor = np.trace(mtensor, axis1=index1, axis2=index2).reshape((newdim, newdim))
    return mtensor


#--------------------------------------------------------------------------

if __name__ == '__main__':
    rho = np.kron(np.kron([[0.5, 0.5], [0.5, 0.5]], [
                  [1, 0], [0, 0]]), [[0, 0], [0, 1]])
    print(partialtrace(rho, [1, 2]))
    # with h5py.File(archivepath, 'a') as archive:
    #circuit = QCircExp(archive, '2QPT_SWAP_ibmqx2')
    # rhoplot(circuit.kraus[-1])
    #condell(archive, '1QPT', ['X', 'X'], mateng)
    #densityplot(archive, 'ProcessTensor3_simulator')
    #comparechiplot('1QPT_U', circuits=9, method='diamond')
