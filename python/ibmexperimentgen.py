# -*- coding: utf-8 -*-
# @Author: Helios
# @Date:   2017-08-21 15:31:32
# @Last Modified by:   Helios
# @Last Modified time: 2017-08-22 17:52:39

import os
import sys
import h5py
import numpy as np
from time import gmtime
import ibmtomography as itm
import ibmupload as iud

#--------------------------------------------------------------------------
# COMPLETED
#--------------------------------------------------------------------------


#--------------------------------------------------------------------------
# TO BE IMPLEMENTED
#--------------------------------------------------------------------------



#--------------------------------------------------------------------------
# CURRENTLY WORKING ON
#--------------------------------------------------------------------------


# generates and saves the base code to run a given number of circuits with 
# a given number of random unitaries. Save file is as 'basename'+'_U' + 'units' + '_' + 'circuititer'
def randomUgen(basename, units, circuits=10):
    if os.path.isdir('C:\\Users\\Helios'):
        basepath = 'C:\\Users\\Helios\\Documents\\Projects\\Uni\\2017\\markovpy\\tomogcircuits'
    else:
        basepath = 'C:\\Users\\Joshua\\Documents\\Projects\\Uni\\2017\\Research\\markovpy\\tomogcircuits'

    # retrieve base qasm for all new files
    with open(basepath + '\\base.qasm', 'r') as base:
        baseqasm = base.read()
    # generate random circuits
    for i in range(0,circuits): 
        circuit = baseqasm
        for j in range(0, units):
            randset = np.random.uniform(size=(3,))
            circuit += 'u3({0}*pi, {0}*pi,{0}*pi) q[0]; \n'.format(*randset)
        filename = basepath + '\\' + basename +'_U{}_C{}'.format(units, i+1)  +'.qasm'
        with open(filename, 'w') as new_code:
            new_code.write(circuit)

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
        filename = basepath + '\\' + basename +'_H{}'.format(j+1)  +'.qasm'
    with open(filename, 'w') as new_code:
            new_code.write(circuit)
#--------------------------------------------------------------------------



if __name__ == '__main__':
    # for i in range(1,11):
    #     Hgen('1QPT_Hadamard3', gates=i)
    #     Hgen('1QPT_Hadamard4', gates=i)
    # add newly made files to file queue
    if os.path.isdir('C:\\Users\\Helios'):
        base = 'C:\\Users\\Helios\\Documents\\Projects\\Uni\\2017\\markovpy\\tomogcircuits\\'
    else:
        base = 'C:\\Users\\Joshua\\Documents\\Projects\\Uni\\2017\\Research\\markovpy\\tomogcircuits\\'

    with open('filequeue.txt', 'a') as fq:
        for i in range(1,11):
            line1 = base + '1QPT_Hadamard3_H{}.qasm;ibmqx2;process;0\n'.format(i) 
            #line2 = base + '1QPT_U{0}_C{1}.qasm;simulator;process;0\n'.format(i,j) 
            fq.write(line1)
            #fq.write(line2)