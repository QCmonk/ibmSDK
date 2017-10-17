# -*- coding: utf-8 -*-
# @Author: Helios
# @Date:   2017-08-21 15:31:32
# @Last Modified by:   Helios
# @Last Modified time: 2017-10-17 12:02:06

import os
import sys
import h5py
import numpy as np
from pyswarm import pso
import ibmupload as iud
from time import gmtime
import ibmtomography as itm
from scipy.optimize import minimize


# gate definition dictionary - add custom gates as needed

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
# Sdg gate
_Sdg = np.matrix([[1, 0], [0, -1j]])
# T gate
_T = np.matrix([[1, 0], [0, (1 + 1j)/np.sqrt(2)]])
# Tdg gate
_Tdg = np.matrix([[1, 0], [0, (1 - 1j)/np.sqrt(2)]])
# CNOT gate
_CX = np.matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
# zero composition
_P1 = np.matrix([[1, 0], [0, 0]])
# one composition
_P2 = np.matrix([[0, 0], [0, 1]])

# gate dictionary
gatedict = {'h':   _H,
            'id':  _ID,
            'x':   _X,
            'y':   _Y,
            'z':   _Z,
            't':   _T,
            'tdg': _Tdg,
            's':   _S,
            'sdg': _Sdg,
            'cx':  _CX,
            'p1':  _P1,
            'p2':  _P2}


#--------------------------------------------------------------------------
# COMPLETED
#--------------------------------------------------------------------------


# generates and saves the base code to run a given number of circuits with
# a given number of random unitaries. Save file is as 'basename'+'_U' +
# 'units' + '_' + 'circuititer'
def randomUgen(basename, units, circuits=10):
    if os.path.isdir('C:\\Users\\Helios'):
        basepath = 'C:\\Users\\Helios\\Documents\\Projects\\Uni\\2017\\markovpy\\tomogcircuits'
    else:
        basepath = 'C:\\Users\\Joshua\\Documents\\Projects\\Uni\\2017\\Research\\markovpy\\tomogcircuits'

    # retrieve base qasm for all new files
    with open(basepath + '\\base.qasm', 'r') as base:
        baseqasm = base.read()
    # generate random circuits
    for i in range(0, circuits):
        circuit = baseqasm
        for j in range(0, units):
            randset = np.random.uniform(size=(3,))
            circuit += 'u3({0}*pi, {0}*pi,{0}*pi) q[0]; \n'.format(*randset)
        filename = basepath + '\\' + basename + \
            '_U{}_C{}'.format(units, i+1) + '.qasm'
        with open(filename, 'w') as new_code:
            new_code.write(circuit)

# generates sequential hadamard gate circuits


def Hgen(basename, gates, inputstate='0'):
    if os.path.isdir('C:\\Users\\Helios'):
        basepath = 'C:\\Users\\Helios\\Documents\\Projects\\Uni\\2017\\markovpy\\tomogcircuits'
    else:
        basepath = 'C:\\Users\\Joshua\\Documents\\Projects\\Uni\\2017\\Research\\markovpy\\tomogcircuits'

    # retrieve base qasm for all new files
    with open(basepath + '\\base.qasm', 'r') as base:
        baseqasm = base.read()
    circuit = baseqasm
    for j in range(0, gates):
        circuit += 'h q[0];\n'
        filename = basepath + '\\' + basename + '_H{}'.format(j+1) + '.qasm'
    with open(filename, 'w') as new_code:
        new_code.write(circuit)

# return required pauli spin matrix


def pauli(num):
    if num == 1:
        return gatedict['x']
    elif num == 2:
        return gatedict['y']
    elif num == 3:
        return gatedict['z']
    else:
        return gatedict['id']

# general single qubit operator


def ibmu3(theta, phi, lmbda):
    return np.asarray([[np.exp(-1j*(phi+lmbda)/2)*np.cos(theta/2), -np.exp(-1j*(phi-lmbda)/2)*np.sin(theta/2)],
                       [np.exp(1j*(phi-lmbda)/2)*np.sin(theta/2), np.exp(1j*(phi+lmbda)/2)*np.cos(theta/2)]])



# generates the 1 qubit map spanning set - no obvious way to do this
# cleanly so here we are
def zbasegen():
    # base element
    zbase = [gatedict['id']]
    for i in range(1, 4):
        baseelp = (np.eye(2) + 1j*pauli(i))/np.sqrt(2)
        baseeln = (np.eye(2) - 1j*pauli(i))/np.sqrt(2)
        # append positive and negative components
        zbase.append(baseelp)
        zbase.append(baseeln)

    z8 = (np.eye(2) + pauli(1)*1j/np.sqrt(2) + pauli(2)*1j/np.sqrt(2))/np.sqrt(2)
    z9 = (np.eye(2) + pauli(1)*1j/np.sqrt(2) + pauli(3)*1j/np.sqrt(2))/np.sqrt(2)
    z10 = (np.eye(2) + pauli(2)*1j/np.sqrt(2) + pauli(3)*1j/np.sqrt(2))/np.sqrt(2)
    zbase.append(z8)
    zbase.append(z9)
    zbase.append(z10)

    return zbase


#--------------------------------------------------------------------------
# TO BE IMPLEMENTED
#--------------------------------------------------------------------------


#--------------------------------------------------------------------------
# CURRENTLY WORKING ON
#--------------------------------------------------------------------------

def spansetgen():
    if os.path.isdir('C:\\Users\\Helios'):
        basepath = 'C:\\Users\\Helios\\Documents\\Projects\\Uni\\2017\\markovpy\\tomogcircuits'
    else:
        basepath = 'C:\\Users\\Joshua\\Documents\\Projects\\Uni\\2017\\Research\\markovpy\\tomogcircuits'

    # retrieve base process tensor file
    with open(basepath + '\\ProcessTensorBase.qasm', 'r') as base:
        baseqasm = base.read()
    circuit = baseqasm

    # generate the spanning set for single qubit maps
    zbase = zbasegen()

    ubase = zconvert(zbase)

    # for each operator, determine the 3 parameters required for arbitary u3

    # add the string combination for that to a new file


def zconvert(zbase):
    # initialise unitary basis equivalents
    ubase = []

    for el in zbase:
        xopt, fopt =swarmfind(el)
        if fopt > 0.5:
            print(el*el.conj().T)
        ubase.append(xopt)
    # convert to string gate equivalents
    ustring = []

# define objective funtion to optimise
def unitaryopt(params, *args):
    op = np.reshape(args, [2,2])
    theta = params[0]*np.pi
    phi = params[1]*np.pi
    lmbda = params[2]*np.pi
    val = np.linalg.norm(ibmu3(theta, phi, lmbda) - op, ord=2)
    return val


# uses particle swarm optimisation to find the best 3 parameters
def swarmfind(op):
    # define upper and lower bounds
    lb = [-2,-2,-2]
    ub = [2,2,2]
    xopt, fopt = pso(func=unitaryopt, lb=lb, ub=ub, args=op, swarmsize=400, phip=1.0)
    return xopt, fopt





#--------------------------------------------------------------------------


if __name__ == '__main__':
    spansetgen()
    #swarmfind(gatedict['id'])
    # with open('filequeue.txt', 'a') as fq:
    # for i in range(1,11):
    #     line1 = base + '1QPT_Hadamard3_H{}.qasm;ibmqx2;process;0\n'.format(i)
    #     #line2 = base + '1QPT_U{0}_C{1}.qasm;simulator;process;0\n'.format(i,j)
    #     fq.write(line1)
    #     #fq.write(line2)
