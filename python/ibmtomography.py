# -*- coding: utf-8 -*-
# @Author: Helios
# @Date:   2017-07-13 14:20:04
# @Last Modified by:   Helios
# @Last Modified time: 2017-09-20 18:16:04


import os
import sys
import h5py
import numpy as np
import matlab.engine
import multiprocessing as mp
from itertools import product
from itertools import combinations



# path to experiment data file
if os.path.isdir('C:\\Users\\Helios'):
    archivepath = 'C:\\Users\\Helios\\Documents\\Projects\\Uni\\2017\\markovpy\\archive.hdf5'
else:
    archivepath = 'C:\\Users\\Joshua\\Documents\\Projects\\Uni\\2017\\Research\\markovpy\\archive.hdf5'


#--------------------------------------------------------------------------
# COMPLETED
#--------------------------------------------------------------------------

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


# returns density operator series term constructed out of the Pauli
# operator basis set. Defaults to Matlab compatible output
def pauliexpand(pauli, mtlb=True):

    # initialise gate
    gate = 1.0

    for el in pauli:
        if el == '0':
            gate = np.kron(gate, gatedict['x'])
        elif el == '1':
            gate = np.kron(gate, gatedict['y'])
        elif el == '2':
            gate = np.kron(gate, gatedict['z'])
        elif el == '3':
            gate = np.kron(gate, gatedict['id'])

    if mtlb:
        return matlab.double(gate.tolist(), is_complex=True)
    else:
        return gate

# zeroes ones in binary string that contribute to normalisation condition


def stokezero(binstr, meas):
    relones = iter(index for index, item in enumerate(meas) if item == '3')
    binstr = list(binstr)
    for index in relones:
        binstr[index] = '0'
    return ''.join(binstr)


# calculates stokes parameter given measurement axis for now
def stokes(labels, probs, measureind):
    # get number of qubits
    qbits = len(measureind)
    # create dictionary using label/prob pairs
    results = {key[:-qbits-1:-1]: value for (key, value) in zip(labels, probs)}
    # create all permutations of binary string
    combs = iter(''.join(seq) for seq in product(['0', '1'], repeat=qbits))
    # iterate over all possible permutations
    stokes = 0.0
    for item in combs:
        # account for normalisation modifier
        binstr = stokezero(item, measureind)
        # probability multiplier
        if binstr.count('1') % 2 == 1:
            mult = -1
        else:
            mult = 1
        # attempt to access probability value for outcome using dictionary
        try:
            prob = mult*results[item]
        # key error for zero probability event (will only happen on simulator I
        # imagine)
        except KeyError:
            prob = 0.0
        stokes += prob

    return stokes


# useful and fast function to flatten list of lists
flatten = lambda l: [item for sublist in l for item in sublist]



# zeroes ones in binary string that contribute to normalisation condition
def stokezero(binstr, meas):
    relones = iter(index for index, item in enumerate(meas) if item == '3')
    binstr = list(binstr)
    for index in relones:
        binstr[index] = '0'
    return ''.join(binstr)


# construct identity deviation generator
def idperms(measind):
    # number of measurements to be swapped out
    replace = measind.count('2')
    replaced = 0

    # find locations of all z axis measurements
    zloc = list(index for index, item in enumerate(measind) if item == '2')
    # compute all possible permutations of swap procedure
    chk = [list(combinations(zloc, bits)) for bits in range(1, replace+1)]
    # flatten list for easy indexing
    chk = flatten(chk)
    # turn original into list for item replacement
    measind = list(measind)
    while replaced < len(chk):
        # copy original for manipulation
        yielditem = measind[:]
        for i in chk[replaced]:
            # replace iteration targets
            yielditem[i] = '3'
        # yield the new measurement sring
        yield ''.join(yielditem)
        # increment generator
        replaced += 1


# perform state tomography for a given circuit using stokes parameter
# method (no linear inversion)
def statetomography_ST(archive, circuit):
    statecirc = archive[circuit]
    # retrieve number of qubits and hence dimensionality of system
    qubits = len(statecirc.attrs['qubits'])

    # iterate through each circuit permutation
    meas = []
    stoke = []
    for perm in statecirc:
        # skip raw qasm file
        if perm == 'raw_qasm':
            continue
        # retrieve measurement basis for permutation (order is 0:x 1:y z:2)
        meas.append(perm[-qubits:])
        label = statecirc[perm]['labels']
        probs = statecirc[perm]['values']
        stoke.append(stokes(label, probs, meas[-1]))
        # get identity stokes as well
        if '2' in meas[-1]:
            iddevs = idperms(meas[-1])
            for dev in iddevs:
                meas.append(dev)
                stoke.append(stokes(label, probs, meas[-1]))

    # construct density matrix
    density = np.zeros((2**qubits, 2**qubits), dtype='complex128')

    for index, param in enumerate(stoke):
        density += param*pauliexpand(meas[index])
    return density/(2**qubits)

# compute trace distance of two matrices
def tracedist(m1, m2):
    dist = 0.5*np.trace(np.sqrt(np.asmatrix(m1-m2).H*(m1-m2)))
    return dist


# perform state tomography using maximun likelihood (requires matlab CVX)
def statetomography_ML(archive, circuit, mateng, goal='state'):
    # create container for target circuit
    statecirc = archive[circuit]
    # retrieve number of qubits and hence dimensionality of system
    qubits = len(statecirc.attrs['qubits'])

    if goal == 'state':
        # generate expectation values from dataset
        meas = []
        stoke = []
        for perm in statecirc:
            # skip raw qasm file
            if perm == 'raw_qasm':
                continue
            # retrieve measurement basis for permutation (order is 0:x 1:y z:2)
            meas.append(perm[-qubits:])
            label = statecirc[perm]['labels']
            probs = statecirc[perm]['values']
            stoke.append(stokes(label, probs, meas[-1]))
            # get identity stokes as well
            if '2' in meas[-1]:
                iddevs = idperms(meas[-1])
                for dev in iddevs:
                    meas.append(dev)
                    stoke.append(stokes(label, probs, meas[-1]))

        # map Pauli basis states into measure list and convert to Matlab
        # supported type
        meas = list(map(pauliexpand, meas[:-1]))
        # perform type conversion on stokes parameters
        stoke = list(map(float, stoke[:-1]))
        # perform MLE with Matlab CVX engine
        density = mateng.statetomography(meas, stoke)
        return density

    # perform state tomgraphy on preperation circuit sets
    elif goal == 'process':
        # initialise density matrix set
        density_set = []
        # length of tomography identifier string
        strlen = 1 + 2*qubits
        # preperation state combinations
        prepcombs = product(['0', '1', '2', '3'], repeat=qubits)
        # iterate over preperation states
        for comb in prepcombs:
            # reset variables for each prep permutation
            meas = []
            stoke = []
            permstring = '_' + ''.join(str(t) for t in comb) + '_'
            for perm in statecirc:
                if perm == 'raw_qasm':
                    continue
                # perform regular state tomography on permutation
                if permstring in perm:
                    # retrieve measurement basis for permutation (order is 0:x
                    # 1:y z:2)
                    meas.append(perm[-qubits:])
                    label = statecirc[perm]['labels']
                    probs = statecirc[perm]['values']
                    stoke.append(stokes(label, probs, meas[-1]))
                    # get identity stokes as well
                    if '2' in meas[-1]:
                        iddevs = idperms(meas[-1])
                        for dev in iddevs:
                            meas.append(dev)
                            stoke.append(stokes(label, probs, meas[-1]))
            # map Pauli basis states into measure list and then convert to
            # Matlab supported type
            meas = list(map(pauliexpand, meas[:-1]))
            # perform type conversion on stokes parameters
            stoke = list(map(float, stoke[:-1]))

            # perform MLE with Matlab CVX engine
            density = mateng.statetomography(meas, stoke)
            density_set.append([permstring[:-1], density])

        return density_set


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
    rax.set_zlim3d([0, dzr.max()])
    iax.set_zlim3d([0, dzr.max()])

    inorm = colors.Normalize(dzi.min(), dzi.max())
    icolours = cm.jet(inorm(dzi))

    # plot image
    rax.bar3d(x, y, z, dx, dy, dzr, color=rcolours)
    iax.bar3d(x, y, z, dx, dy, dzi, color=icolours)
    plt.ticklabel_format(style='sci', axis='z', scilimits=(0, 0))
    plt.show()


# compute Kraus operators given process matrix and basis set
def mapcompute(chi, operatorbasis):
    from scipy.linalg import eigh
    # compute diag of process matrix
    D, V = eigh(chi)
    dim = int(np.log2(V.shape[0]))
    # compute unitary
    U = V.conj().T
    # initialise Kraus set
    kraus = np.asarray([np.zeros((dim, dim), dtype='complex128')]*(dim**2))
    # compute Kraus operators
    for i in range(0, dim**2):
        for j in range(0, dim**2):
            kraus[i] += np.sqrt(D[i])*operatorbasis[j]*U[i][j]
    return kraus




# converts the kraus represention of a channel to the Choi form
def kraus2choi(kraus, rep='choi', targets=1):
    # compute dimension of system
    if targets==1:
        dim = kraus[0].shape[0]
    else:
        dim = (2**targets)
    # initialise A form
    loui = np.zeros([dim**2]*2, dtype='complex128')
    # compute liouville representation of Kraus operator
    for i in kraus:
        if targets == 1:
            loui += np.kron(i.conj(), i)
        else:
            loui += np.kron(np.kron(i, np.eye(2)).conj(),np.kron(i, np.eye(2)))
    # return appropriate representation
    if rep == 'loui':
        return loui
    else:
        # permute to Choi representation
        choi = loui.reshape([dim, dim, dim, dim])
        choi = np.transpose(choi, [0, 2, 1, 3])
        choi = choi.reshape([dim**2, dim**2])
        return choi



# compute effect of a CPTP map
def CPTP(kraus, rho):
    
    nrho = np.zeros(rho.shape, dtype='complex128')
    for i in kraus:
        nrho += np.dot(i, np.dot(rho, i.conj()))
    return nrho




# formats a numpy array into a matlab compatible vector structure
def numpy2matlabvec(array):
    from functools import partial
    # define mappable function with complex argument
    mapmatlab = partial(matlab.double, is_complex=True)
    # perform conversion
    return list(map(mapmatlab, [item.tolist() for item in array]))


#
def numpyarr2matlab(array):
    array = np.asarray(array)
    # perform conversion
    return matlab.double(array.tolist(), is_complex=True)


# generates tensor products of arrayset according to combination
def kronjob(arrayset, combination):
    # initialise array
    matrix = 1
    # construct appropriate tensor product
    for el in combination:
        matrix = np.kron(matrix, arrayset[int(el)])
    return matrix



# generates preperation basis states and converts to matlab compatible
# type if requested

def prepbasis(qnum, mtlb=False):
    # define base set
    prepset = np.asarray([[[1, 0], [0, 0]],
                          [[0, 0], [0, 1]],
                          [[1/2, 1/2], [1/2, 1/2]],
                          [[1/2, -1j/2], [1j/2, 1/2]]])

    # determine all combinations
    combs = iter(''.join(seq)
                 for seq in product(['0', '1', '2', '3'], repeat=qnum))
    densitybasis = []
    # construct density basis
    for item in combs:
        densitybasis.append(kronjob(prepset, item))
    densitybasis = np.asarray(densitybasis)

    # convert to matlab compatible type if required
    if mtlb:
        return numpy2matlabvec(densitybasis)
    else:
        return densitybasis



# generates an operator spanning set for any system with qnum qubits
def opbasis(qnum, mtlb=False):
    # define operator basis set (Pauli set in this case)
    opset = np.asarray([gatedict['id'],
                        gatedict['x'],
                        gatedict['y'],
                        gatedict['z']])

    # determine all combinations
    combs = iter(''.join(seq)
                 for seq in product(['0', '1', '2', '3'], repeat=qnum))
    operbasis = []
    # construct density basis
    for item in combs:
        operbasis.append(kronjob(opset, item))
    operbasis = np.asarray(operbasis)

    # convert to matlab if required
    if mtlb:
        return numpy2matlabvec(operbasis)
    else:
        return operbasis


# define beta worker
def betawork(densitybasis, operatorbasis, partition):
    import matlab.engine
    mateng = matlab.engine.start_matlab()
    betaseg = mateng.betaopt(densitybasis, operatorbasis, partition)
    return np.asarray(betaseg)


# determines division of labor between cores
def labourdiv(dim, cores):
    # compute full workload range (+1 because of matlab indexing)
    wload = list(range(1, dim**2+1))
    # compute ~number of iterations per worker
    witer = len(wload)//cores
    # yield cores ranges until wload is empty
    for i in range(0, cores):
        # export all remaining work
        if i == cores-1:
            yield matlab.double([min(wload[i*witer:]), max(wload[i*witer:])])
        else:
            yield matlab.double([min(wload[i*witer:(i+1)*witer]), max(wload[i*witer:(i+1)*witer])])


# perform process tomography for a given circuit, assuming all circuits
# are present
def processtomography(archive, circuit, mateng, parallel=True):
    # get number of qubits to construct basis for
    qnum = len(archive[circuit].attrs['qubits'][()])
    dim = int(2**qnum)
    # define standard basis set for CPTP map reconstruction
    # (|0>,|1>,|+>,|+i>), only single qubit at the moment
    densitybasis = prepbasis(qnum, mtlb=True)
    # define Kraus operator basis set
    operatorbasis = opbasis(qnum, mtlb=True)
    # retrieve density operators after map has acted on spanning set
    reconstruct = []
    tomogpath = archive[circuit]['tomography_ML']
    for state in tomogpath:
        reconstruct.append(matlab.double(
            tomogpath[state][()].tolist(), is_complex=True))

    # get number of available cores:
    cores = int(os.environ['NUMBER_OF_PROCESSORS'])
    # parallelise problem if requested and can be done in system
    betapath = 'Support_Group/Process_beta{}'.format(qnum)
    if parallel and cores > 1 and qnum > 1:
        # It is less clear how to use Matabs inbuilt parallelisation for this problem
        # so I opted for implementing it in python. Modify as you see fit.
        
        # perform labmda optimisation sequentially
        lam = mateng.lambdaopt(reconstruct, densitybasis)

        # retrieve beta matrix if already computed else re-compute
        if betapath in archive:
            beta = archive['Support_Group']['Process_beta{}'.format(int(qnum))][()]
        else:
            print('------- DISTRIBUTING HYPERDIMENSIONAL OPTIMISATION PROBLEM ACROSS {} CORES -------'.format(cores))
            # perform beta optimisation in parallel if required
            # determine work division
            wrange = labourdiv(dim, cores)
            # worker list
            pool = mp.Pool(processes=cores)
            workers = []
            for i, partition in enumerate(wrange):
                workers.append(pool.apply_async(betawork, args=(
                    densitybasis, operatorbasis, partition)))
            pool.close()
            pool.join()

            # retrieve results
            beta = np.zeros((dim**4, dim**4), dtype='complex128')
            for employee in workers:
                beta += employee.get()

        # compute chi matrix sequentially
        chi = np.asarray(mateng.chiopt(
            numpyarr2matlab(beta), operatorbasis, lam))

    else:
        
        chi = np.asarray(mateng.processtomography(
            reconstruct, densitybasis, operatorbasis))
    # output result
    return chi, np.asarray(operatorbasis), np.asarray(densitybasis)



# computes the beta matrices for process tomography up to the maximum number of qubits 
# and stores in support groups - GREATLY reduces compute time for process tomography
def betacompute(qubits=2, parallel=True):
    mateng = matlab.engine.start_matlab()
    betabase = 'Process_beta'
    # get number of available cores:
    cores = int(os.environ['NUMBER_OF_PROCESSORS'])
    with h5py.File(archivepath, 'a') as archive:
        try:
            # check if support group exists
            if 'Support_Group' not in archive:
                archive.create_group('Support_Group')
            for i in range(1, qubits+1):

                # define standard basis set for CPTP map reconstruction
                # (|0>,|1>,|+>,|+i>), only single qubit at the moment
                densitybasis = prepbasis(i, mtlb=True)
                # define Kraus operator basis set
                operatorbasis = opbasis(i, mtlb=True)
                # compute the dimension of the problem
                dim = 2**i
                # compute beta matrix
                betaname = betabase + '{}'.format(i)
                # delete beta matrix if it already exists
                if betaname in archive['Support_Group']:
                    archive['Support_Group'].__delitem__(betaname)

                if parallel and i>1:
                    print('------- DISTRIBUTING HYPERDIMENSIONAL OPTIMISATION PROBLEM ACROSS {} CORES -------'.format(cores))
                    # perform beta optimisation in parallel
                    # determine work division
                    wrange = labourdiv(dim, cores)
                    # worker list
                    pool = mp.Pool(processes=cores)
                    workers = []
                    for i, partition in enumerate(wrange):
                        workers.append(pool.apply_async(betawork, args=(
                            densitybasis, operatorbasis, partition)))
                    pool.close()
                    pool.join()
                    # retrieve results
                    beta = np.zeros((dim**4, dim**4), dtype='complex128')
                    for employee in workers:
                        beta += employee.get() 
                else:
                    partition = next(labourdiv(dim, 1))
                    beta = mateng.betaopt(densitybasis, operatorbasis, partition)

                # add computed matrix to archive
                archive['Support_Group'].create_dataset(betaname, data=np.asarray(beta))

        except Exception as e:
            print(e)
            archive.flush()
            archive.close()
            exit()


# computes the partial trace of m \in \mathcal{H^n}, tracing out subsystems in sys
# why must you make my life so difficult numpy?
def partialtrace(m, sys):
    # type enforcement
    m = np.asarray(m)
    # sort subsystems
    sys = sorted(sys)
    # get tensor dimensions
    qnum = int(round(np.log2(len(m))))
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




# computes the output of a quantum channel given the input rho and the maps B form
def choimap(choi, rho):
    # compute dimensionality of system
    dimr = len(rho)
    dimc = len(choi)
    # compute joint input state
    hinput = np.kron(np.eye(dimr), rho.T)
    # compute subsystems to trace out
    subsys = [i for i in range(dimr//2, dimr)]
    # trace out input dimensions for final state
    return partialtrace(np.dot(hinput, choi), subsys)



# DANGER FUNCTIONALITY, USE AT OWN RISK
# deletes a particular item (e.g Data_Group) from every group in root directory
def _blockdelete(group):
    with h5py.File(archivepath, 'a') as archive: 
        for item in archive:
            if group in archive[item]:
                archive[item].__delitem__(group)

# add default shots if parameter does not exist (assume lowest shot number ever used)
def _shotadd(shotnum=4096):
    with h5py.File(archivepath, 'a') as archive:
        for item in archive:
            if archive[item].attrs.__contains__('shots'):
                continue
            else:
                archive[item].attrs['shots'] = shotnum



# computes distance of a given process tensor from that of the closest markovian process
# using the chosen distance metric: quantum relative entropy or tracenorm (only included for thoroughness)
def markovdist(ptensor, subsystems, metric='qre'):
    # compute dimensionality of sysem (corresponds to number of timesteps)
    mptensor = 1.0
    for item in subsystems:
        mptensor = np.kron(mptensor,partialtrace(ptensor, item))
    return mptensor
    # compute distance using chosen metric
    if metric == 'qre':
        # compute quantum relative entrop
        from scipy.linalg import logm
        return np.trace(np.dot(ptensor,(logm(ptensor) - logm(mptensor))))
    elif metric == 'trnm':
        # compute trace norm
        return tracenorm(ptensor, mptensor)
    else:
        print('Unknown metric')
        return np.inf



#--------------------------------------------------------------------------
# TO BE IMPLEMENTED
#--------------------------------------------------------------------------


#--------------------------------------------------------------------------
# CURRENTLY WORKING ON
#--------------------------------------------------------------------------


#--------------------------------------------------------------------------

if __name__ == '__main__':
    #_blockdelete('Data_Group')
    testing = np.kron([[0.5,-0.5],[-0.5,0.5]],[[0,0],[0,1]])
    testing = np.kron(testing, testing)
    print( markovdist(testing, [[0,1], [2,3]]) - testing)
    # with h5py.File(archivepath, 'a') as archive:
    # #     # archive.__delitem__('hadamardq0_simulator')
    #      try:
    #         import ibmanalysis as ibma
    #         rho = np.kron(np.asarray([[1,0],[0,0]]),np.asarray([[0,0],[0,1]]))
    #         choi = archive['2QPT_CX_ibmqx2']['Data_Group']['Choi_matrix']
    #         ibma.rhoplot(rho)
    #         ibma.rhoplot(choimap(choi, rho))
    #      except Exception as e:
    #         print(e)
    #         archive.flush()
    #         archive.close()
    #         exit()

