# -*- coding: utf-8 -*-
# @Author: Helios
# @Date:   2017-10-10 18:41:02
# @Last Modified by:   Helios
# @Last Modified time: 2017-10-17 11:59:32

import os
import sys
import h5py
import matplotlib
import numpy as np
import matlab.engine
from time import gmtime
import ibmtomography as itm
import ibmanalysis as ianal
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



# identity matrix
_ID = np.matrix([[1, 0], [0, 1]])
# X gate
_X = np.matrix([[0, 1], [1, 0]])
# Z gate
_Z = np.matrix([[1, 0], [0, -1]])
# Hadamard gate
_H = (1/np.sqrt(2))*np.matrix([[1, 1], [1, -1]])
# Y Gate
_Y = np.matrix([[0, -1j], [1j, 0]])
# S gate
_S = np.matrix([[1, 0], [0, 1j]])
# T gate
_T = np.matrix([[1, 0], [0, (1 + 1j)/np.sqrt(2)]])
# CNOT gate
_CX = np.matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
# random gate
_FF = ibmu3(0.22*np.pi, 0.65*np.pi, 0.3*np.pi)
# gate dictionary
gatedict = {'H':   _H,
            'I':   _ID,
            'X':   _X,
            'Y':   _Y,
            'Z':   _Z,
            'T':   _T,
            'S':   _S,
            'CX':  _CX,
            'FF':  _FF}


def ibmu3(theta, phi, lmbda):
    return np.asarray([[np.exp(-1j*(phi+lmbda)/2)*np.cos(theta/2), -np.exp(-1j*(phi-lmbda)/2)*np.sin(theta/2)],
                       [np.exp(1j*(phi-lmbda)/2)*np.sin(theta/2), np.exp(1j*(phi+lmbda)/2)*np.cos(theta/2)]])



#--------------------------------------------------------------------------
# CURRENTLY WORKING ON
#--------------------------------------------------------------------------


# quantum state class, for simulating quantum systems
class QCirc(object):

    def __init__(self, qubits, init='default'):
        # number of qubits in system
        self.qubits = qubits
        # initial state of individual qubits
        self.init = init
        # store a matlab engine instance
        self.mateng = matlab.engine.start_matlab()
        # initiate
        self.circinit()

    # state initialisation
    def circinit(self):

        # define base states of inidiviudal qubits, defaults to input
        if self.init == 'default':
            basestate = np.asarray([[1, 0], [0, 0]])
        elif self.init == 'super':
            basestate = np.asarray([[1, 1], [1, 1]])/2
        else:
            basestate = self.init

        # initialise quantum state with required number of qubits using density
        # operator formalism
        state = 1.0
        for i in range(0, self.qubits):
            state = np.kron(state, basestate)

        # store as current state
        self.cstate = state

    # applies choi map to state targets
    def choiapply(self, choi, target):
        # perform type conversion
        choimat   = itm.numpyarr2matlab(choi)
        statemat  = itm.numpyarr2matlab(self.cstate)
        targetmat = [i+1 for i in target]
        # compute map with change to qubit indexing
        self.cstate = np.asarray(self.mateng.mapcompute(choimat, statemat, targetmat))

    def markov(self, qubits):
        #compute markov distance
        self.mdist, self.markov = ianal.markovdist(self.cstate, qubits, gen=True)

    def nmarkov(self, qre):
        # perform type conversion
        markovmat = itm.numpyarr2matlab(self.markov)
        statemat  = itm.numpyarr2matlab(self.cstate)

        # compute choi map that corrects nonMarkovianity 
        self.corrchoi = np.asarray(self.mateng.correlationmap(markovmat, statemat, qre))




# compute the process tensor given the two time steps

def statecompute(archive, g1, g2, plot=False, save=False, noise='sim'):
    # extract choi states of gates
    if noise=='sim':
        cnot = archive['2QPT_CX_simulator']['Data_Group']['Choi_matrix'][()]
        had = archive['1QPT_H_simulator']['Data_Group']['Choi_matrix'][()]
        gate1 = archive['1QPT_' + g1 + '_simulator']['Data_Group']['Choi_matrix'][()]
        gate2 = archive['1QPT_' + g2 + '_simulator']['Data_Group']['Choi_matrix'][()]

    elif noise=='actual':
        cnot = archive['2QPT_CX_ibmqx4']['Data_Group']['Choi_matrix'][()]
        had = archive['1QPT_H_ibmqx4']['Data_Group']['Choi_matrix'][()]
        gate1 = archive['1QPT_' + g1 + '_ibmqx4']['Data_Group']['Choi_matrix'][()]
        gate2 = archive['1QPT_' + g2 + '_ibmqx4']['Data_Group']['Choi_matrix'][()]
    elif noise=='none':
        # generate choi states of perfect maps
        cnot = itm.kraus2choi(gatedict['CX'])
        had = itm.kraus2choi(gatedict['H'])
        gate1 = itm.kraus2choi(gatedict[g1])
        gate2 = itm.kraus2choi(gatedict[g2])

    # initialise circuit class
    approxcirc = QCirc(qubits=4, init='default') 

    # generate entanled states
    approxcirc.choiapply(had, [1])
    approxcirc.choiapply(cnot, [0,1])
    approxcirc.choiapply(had, [3])
    approxcirc.choiapply(cnot, [2,3])

    # apply first gate
    approxcirc.choiapply(gate1, [0])

    # apply swap operation
    approxcirc.choiapply(cnot,[0,2])
    approxcirc.choiapply(had, [0])
    approxcirc.choiapply(had, [2])
    approxcirc.choiapply(cnot,[0,2])
    approxcirc.choiapply(had, [2])
    approxcirc.choiapply(had, [0])
    approxcirc.choiapply(cnot,[0,2])

    # apply second gate
    approxcirc.choiapply(gate2, [0])

    # compute markov distance
    approxcirc.markov([[1,2],[0,3]])

    # compute circuit approximation with given non markovianity
    approxcirc.nmarkov(0.22)
    # save correlation matrix if requested
    if save:
        # define corr matrix name
        strname = g1+g2
        # delete choi state if already exists
        if archive['Support_Group'].__contains__('Correlation'+strname):
            archive['Support_Group'].__delitem__('Correlation'+strname)
        archive['Support_Group'].create_dataset('Correlation'+strname, data=approxcirc.corrchoi)
    # plot final state if requested
    if plot:
        ianal.rhoplot(approxcirc.cstate)
        plt.matshow(np.real(approxcirc.corrchoi))
        plt.colorbar()
        plt.show()


def tensortomography(archive, gateset, steps=2, noise='none'):
    # define base string 
    choistate = []
    # compute permuations of base set and extract choi states of base set
    combs = list(product(gateset, repeat=steps))

    # use IBM simulated gates
    if noise=='sim':
        for el in combs:
            basename = 'ProcessTensor'+''.join(el)+'_simulator'
            choistate.append(matlab.double(archive[basename]['tomography_ML'][()].tolist(), is_complex=True))
    elif noise=='none':
    # use perfect gates
        for el in combs:
            choistate.append(matlab.double(pprocessgen2(el[0], el[1]).tolist(), is_complex=True))
    # use actual gates
    elif noise=='actual':
        pass
    else:
        print('------- UNRECOGNISED NOISE PARAMETER: {} -------'.format(noise))

    # retrieve target choi state
    #ptensor = itm.numpyarr2matlab(archive['ProcessTensorUV2_simulator']['tomography_ML'][()])
    ptensor = itm.numpyarr2matlab(pprocessgen2('FF', 'H'))

    # spool up matlab engine
    mateng = matlab.engine.start_matlab()

    # compute lambda set
    lambstoslaughter = np.asarray(mateng.markovmix(choistate, ptensor)[0])
    print(lambstoslaughter)
    print(np.sum(lambstoslaughter))
    # compute approximated process tensor 
    approxtensor = np.zeros([4**steps]*2, dtype='complex128')
    for num,item in enumerate(lambstoslaughter):
        basechoi = np.asarray(choistate[num])
        approxtensor += item*basechoi

    ianal.rhoplot(approxtensor)
    ianal.rhoplot(np.asarray(ptensor))


# generate 2 step process tensor with given parameters and no noise
def pprocessgen2(g1, g2):

    # generate choi states of perfect maps
    cnot = itm.kraus2choi([gatedict['CX']])
    had = itm.kraus2choi([gatedict['H']])
    gate1 = itm.kraus2choi([gatedict[g1]])
    gate2 = itm.kraus2choi([gatedict[g2]])

    # initiate quantum circuit
    ptensor = QCirc(qubits=4, init='default')
    # generate entangled states
    ptensor.choiapply(had, [0])
    ptensor.choiapply(cnot, [0,1])
    ptensor.choiapply(had, [2])
    ptensor.choiapply(cnot, [2,3])

    # apply first gate
    ptensor.choiapply(gate1, [0])

    # apply swap operation
    ptensor.choiapply(cnot,[0,2])
    ptensor.choiapply(had, [0])
    ptensor.choiapply(had, [2])
    ptensor.choiapply(cnot,[0,2])
    ptensor.choiapply(had, [2])
    ptensor.choiapply(had, [0])
    ptensor.choiapply(cnot,[0,2])

    # apply second gate
    ptensor.choiapply(gate2, [0])
    return ptensor.cstate


# correlation map manipulation
def correlationmap(archive, base):
    choicorr = archive['Support_Group']['Correlation'+base][()]
    plt.matshow(np.real(choicorr))
    plt.colorbar()
    plt.show()



#--------------------------------------------------------------------------




if __name__ == '__main__':
    #with h5py.File(archivepath, 'a') as archive:
        #ianal.rhoplot(archive['1QPT_Y_simulator']['Data_Group']['Choi_matrix'][()])
    tensortomography(1,  ['X','Y','Z','I','S','T','H'])
    #statecompute(1, 'X', 'X', plot=False, noise='none')
    #statecompute(archive, 'X', 'X', plot=False, save=True)
    #statecompute(archive, 'X', 'Y', plot=False, save=True)
    #statecompute(archive, 'X', 'Z', plot=False, save=True)
    #statecompute(archive, 'X', 'I', plot=False, save=True)
    #correlationmap(archive, 'XX')
