# -*- coding: utf-8 -*-
# @Author: Helios
# @Date:   2017-07-13 14:20:04
# @Last Modified by:   Helios
# @Last Modified time: 2017-09-01 17:46:11

import os
import sys
import h5py
import numpy as np
from time import gmtime
import ibmtomography as itm
import matplotlib
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
    def __init__(self, target, circuit):
        # store circuit name
        self.name = circuit
        # retrieve and store circuit attributes
        self.device = target.attrs['device'].decode('ascii')
        self.tomography = target.attrs['tomography'].decode('ascii')
        self.qubits = target.attrs['qubits'][()]
        self.complete = target.attrs['complete'][()]
        self.total = target.attrs['total'][()]
        self.rawqasm = target['raw_qasm'][()]

        # check for circuit reconstruction data structures
        if self.complete == self.total:
            # store key information from either
            if self.tomography == 'state':
                self.finalstate = target['tomography_ML'][()]
            elif self.tomography == 'process':
                # extract Chi matrix
                self.chi = np.asarray(target['Data_Group']['Process_matrix'][()])
                # Kraus operators
                self.kraus = np.asarray(target['Data_Group']['Kraus_set'][()])
                # operator basis used in reconstruction
                self.opbasis = np.asarray(target['Data_Group']['Operator_basis_set'][()])
                # preperation basis used in reconstruction
                self.prepbasis = np.asarray(target['Data_Group']['Preperation_basis_set'][()])
                # Choi matrix
                self.choi = np.asarray(target['Data_Group']['Choi_matrix'][()])
                # extract tomography dataset
                self.tomographyset = []
                for state in target['tomography_ML']:
                    self.tomographyset.append(target['tomography_ML'][state][()])



# imports target data from archive and converts to analysis class
def circuitimport(circuit):
    with h5py.File(archivepath, 'a') as archive:
        # create shortcut variable and initialise class
        return QCircExp(archive[circuit], circuit)

#compute trace distance of two matrices using good ol trace norm
def tracedistcrude(m1,m2):
    dist = 0.5*np.trace(np.sqrt(np.asmatrix(m1-m2).H*(m1-m2)))
    return np.real(dist)



# compares actual and simulator results process tomography using some base name
def comparechiplot(basename, circuits=10, comparison='trace'):
    # iterate over different versions and add circuit information to primary data container
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
    if comparison == 'trace':
        for pair in circcollect:
            channeldist = np.append(channeldist,tracedistcrude(pair[0].chi, pair[1].chi))
    elif comparison == 'diamond':
        # import and initialise a matlab engine instance
        import matlab.engine
        mateng = matlab.engine.start_matlab()
        # compute diamond norm for 'equivalent' quantum channels
        for pair in circcollect:
            diffchoi = pair[0].choi - pair[1].choi
            diffchoi = matlab.double(diffchoi.tolist(), is_complex=True)
            channeldist = np.append(channeldist, mateng.dnorm(diffchoi))
    else:
        print('------- UNKNOWN COMPARISON METHOD {}: ABORTING -------'.format(comparison))
        exit()
    gatenum = list(range(1,len(channeldist)+1))
    plt.plot(channeldist)
    plt.title('Channel distance as a function of gate number for circuit: ' + basename)
    plt.grid(True)
    plt.xlabel('Number of applied random unitaries')
    plt.ylabel('Channel distance using ' + comparison + ' norm')
    plt.show()


def tracenormcompare(basename, circuits=10, method='sim'):
    # iterate over different versions and add circuit information to primary data container
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
    if method=='sim':
        for pair in circcollect:
            #print(tracenorm(pair[0].tomographyset[0] - pair[0].tomographyset[1]))
            tracedist = np.append(tracedist, tracenorm(pair[0].tomographyset[0] - pair[0].tomographyset[1]))
    else:
        for pair in circcollect:
            tracedist = np.append(tracedist, tracenorm(pair[1].tomographyset[0] - pair[1].tomographyset[1]))

    gatenum = list(range(1,len(tracedist)+1))
    plt.plot(gatenum, tracedist)
    plt.title('Trace norm as a function of hadamard gates applied for circuit: ' + basename)
    plt.grid(True)
    plt.ylim([1.8,2.1])
    plt.xlabel('Number of applied gates')
    plt.ylabel('distance using tracenorm')
    plt.legend('Actual', 'Ideal')
    plt.show()


def tracenormcomparetmp(circuits=10, method='sim'):
    # iterate over different versions and add circuit information to primary data container
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
        average0 = group[0].tomographyset[0] + group[1].tomographyset[0] + group[2].tomographyset[0] + group[3].tomographyset[0]
        average0 = average0/4.0
        average1 = group[0].tomographyset[1] + group[1].tomographyset[1] + group[2].tomographyset[1] + group[3].tomographyset[1]
        average1 = average1/4.0

        tracedist = np.append(tracedist, tracenorm(average0 - average1))

    gatenum = list(range(1,len(tracedist)+1))
    plt.plot(gatenum, tracedist, label='Actual', linewidth=3.0)
    plt.title('Trace norm as a function of hadamard gates applied for circuit: 1QPT_Hadamard')
    plt.grid(True)
    plt.ylim([1.8,2.05])
    plt.plot( gatenum, [2.0]*len(tracedist), label='Ideal', linewidth=3.0, ls='dashed')
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
    fig.canvas.set_window_title('Process matrix for circuit {}'.format(circuit))
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



#--------------------------------------------------------------------------
# CURRENTLY WORKING ON
#--------------------------------------------------------------------------



#--------------------------------------------------------------------------

if __name__ == '__main__':
    with h5py.File(archivepath, 'a') as archive:
        #chiplot(archive, '2QPT_SWAP_simulator')
        #densityplot(archive, 'ProcessTensor_ibmqx2')
        rho = np.asarray([[1,0],[0,0]])
        circuit = QCircExp(archive['1QPT_H_simulator'], 'Pauli Gate')
        rhoplot(circuit.choi)
