# -*- coding: utf-8 -*-
# @Author: Helios
# @Date:   2017-07-13 14:20:04
# @Last Modified by:   Helios
# @Last Modified time: 2017-08-31 14:08:05


import os
import sys
import h5py
import numpy as np
from time import gmtime
import ibmtomography as itm
import matlab.engine

# python 2 and 3 import of the IBM experience module
if sys.version_info.major > 2:
    from IBMQuantumExperience.IBMQuantumExperience import IBMQuantumExperience
else:
    from IBMQuantumExperience import IBMQuantumExperience


# define user token
token = 'a42ac7aef6ce008f2e7f437a708e64e6683f6dcb0fe043474663b8f4909a40e3357266a496c7c3659c69dc984ba56476d1c6b4c84c4af2cb2d0f7f6d52874cab'

# define config URL explicitly
config = {"url": 'https://quantumexperience.ng.bluemix.net/api'}

# path to experiment data file
if os.path.isdir('C:\\Users\\Helios'):
    archivepath = 'C:\\Users\\Helios\\Documents\\Projects\\Uni\\2017\\markovpy\\archive.hdf5'
else:
    archivepath = 'C:\\Users\\Joshua\\Documents\\Projects\\Uni\\2017\\Research\\markovpy\\archive.hdf5'


#--------------------------------------------------------------------------
# COMPLETED
#--------------------------------------------------------------------------


class TCircGen():
    """
    Implements an iterator that returns the next circuit suffix for a tomography run
    """

    def __init__(self, targets, measureset=[]):
        # store qubits to be reconstructed
        self.targets = targets

        # define the measurement basis set
        if len(measureset) == 0:
            self.measureset = np.asarray(
                ['h q[{0}];\nmeasure q[{0}] -> c[{0}];', 's q[{0}];\nh q[{0}];\nmeasure q[{0}] -> c[{0}];', 'measure q[{0}] -> c[{0}];'])
        else:
            self.measureset = measureset

        # define all combinations of measure set for given number of qubits
        from itertools import product
        self.combs = product(self.measureset, repeat=len(self.targets))
        self.measindex = product(
            list(range(0, len(self.measureset))), repeat=len(self.targets))

    def __iter__(self):
        # initialise circuit iterator number
        self.circnum = 0
        return self

    def __next__(self):

        if self.circnum < len(self.measureset)**len(self.targets):
            # retrieve next circuit configuration
            combination = next(self.combs)
            circuitsuffix = ''

            # generate circuit suffix
            for num, qtarget in enumerate(self.targets):
                circuitsuffix += '\n' + combination[num].format(qtarget)

            # compute measurement indexes applied in current circuit
            circconfig = (len(self.targets)*'{}').format(*next(self.measindex))

            # increment iterator
            self.circnum += 1
            return circuitsuffix, circconfig

        else:
            raise StopIteration


# format target file into appropriate string
def fileparse(filename, filepath=None):
    # retrieve file contents, stripping white space and semicolon delimiter
    if filepath == None:
        with open(os.getcwd() + '\\' + filename, 'r') as f:
            qasm = [m.lstrip() for m in f]
    else:
        with open(filepath + '\\' + filename, 'r') as f:
            qasm = [m.lstrip() for m in f]

    # get rid of all blank lines and return concatenated file
    return ''.join(list(filter(None, qasm)))


# test function for archive
def archivetest(archive):
    for i in archive:
        print(i)


# distributes the execution codes for an uploaded file. Assumes IBM return
# an ordered dictionary D:
def codedist(archive, codes):
    # iterate through archive
    for circuit in archive:
        # check for complete case
        if archive[circuit].attrs['complete'] < archive[circuit].attrs['total']:
            # iterate over first len(codes) circuits with ready status
            for perm in archive[circuit]:
                if perm == 'raw_qasm' or perm == 'tomography_ML' or perm == 'Data_Group':
                    continue

                # check for remaining return codes
                elif len(codes) == 0:
                    break
                # update circuit permutation with return code
                elif archive[circuit][perm].attrs['status'] == 0:
                    archive[circuit][perm].attrs['executionId'] = np.string_(
                        codes[0]['executionId'])
                    archive[circuit][perm].attrs['status'] = 1
                    # remove return code for perm circuit
                    del codes[0]


# adds new circuit to archive
def archiveadd(archive, path, fileparams, shots):
    filepath = os.path.dirname(path)
    filename = os.path.basename(path)
    print('------- ADDING CIRCUIT: {} TO ARCHIVE -------'.format(filename))

    if filename[-5:] != '.qasm':
        print('------- UNRECOGNISED INPUT FILE - ONLY FILES WITH \'.qasm\' EXTENSION ACCEPTED: ABORTING -------')
        # close archive gracefully
        archive.flush()
        archive.close()
        exit()

    # compute uploadable qasm code
    qasm = fileparse(filename, filepath)

    # attempt to add circuit to archive
    try:
        # check if file already in archive (run increase shot number routine)
        dataname = filename[:-5] + '_' + fileparams[0]
        if dataname in archive:
            # check if still waiting on previous results
            if archive[dataname].attrs['complete'] < archive[dataname].attrs['total']:
                print('------- PREVIOUS ITERATION OF CIRCUIT {} IS ONGOING: IGNORING RESUBMISSION COMMAND -------'.format(dataname))
            else:
                print('------- CIRCUIT {} HAS ALREADY BEEN UPLOADED: REPEATING EXPERIMENT/ANALYSIS -------'.format(dataname))
                
                # delete tomography related data
                try:
                    if archive[dataname].attrs['tomography'].decode('ascii') == 'state':
                        archive[dataname].__delitem__('tomography_ML')
                    elif  archive[dataname].attrs['tomography'].decode('ascii') == 'process':
                        archive[dataname].__delitem__('tomography_ML')
                        archive[dataname].__delitem__('Data_Group')
                # fine for now, change to specific error later
                except Exception:
                    pass
                # adjust shot number accordingly
                archive[dataname].attrs['shots'] = archive[dataname].attrs['shots'] + shots 

                # set status for circuit reupload
                for perm in archive[dataname]:
                    # ignore qasm file
                    if perm == 'raw_qasm':
                        continue
                        # set status code for re-upload
                    archive[dataname][perm].attrs['status'] = 0
                    archive[dataname][perm].attrs['executionId'] = 'Null'
                archive[dataname].attrs['complete'] = 0

        # add state tomography circuit
        elif fileparams[1] == 'state':
            newcircuit = archive.create_group(dataname)

            # add appropriate attributes
            # device type to run circuits on
            newcircuit.attrs['device'] = np.string_(fileparams[0])
            # type of tomography that is being performed (state or process)
            newcircuit.attrs['tomography'] = np.string_(fileparams[1])
            # target qubits to measure
            newcircuit.attrs['qubits'] = fileparams[2]
            # number of completed circuits (results retrieved)
            newcircuit.attrs['complete'] = 0
            # total number of circuits to be run
            newcircuit.attrs['total'] = 3**len(fileparams[2])
            # set number of shots
            newcircuit.attrs['shots'] = shots
            # add raw qasm code to group
            newcircuit.create_dataset('raw_qasm', data=np.string_(qasm))

            # generate circuits to be computed
            for circuit, meas in TCircGen(fileparams[2]):
                # create data groups with appropriate format
                grpname = filename[
                    :-5] + '_{}:{}:{}:{}:{}:{}_'.format(*gmtime()[0:6]) + meas
                tomogcirc = newcircuit.create_group(grpname)
                # add status attribute: 0 for ready, 1 for uploaded, 2 for
                # complete
                tomogcirc.attrs['status'] = int(0)
                # add idExecution code (null for now)
                tomogcirc.attrs['executionId'] = np.string_('Null')
                # add specific qasm code
                tomogcirc.create_dataset(
                    'qasm', data=np.string_(qasm + circuit))

        # generate process tomography circuits
        elif fileparams[1] == 'process':
            newcircuit = archive.create_group(dataname)

            # add appropriate attributes
            # device type to run circuits on
            newcircuit.attrs['device'] = np.string_(fileparams[0])
            # type of tomography that is being performed (state or process)
            newcircuit.attrs['tomography'] = np.string_(fileparams[1])
            # target qubits to measure
            newcircuit.attrs['qubits'] = fileparams[2]
            # number of completed circuits (results retrieved)
            newcircuit.attrs['complete'] = 0
            # total number of circuits to be run
            newcircuit.attrs['total'] = (
                3**len(fileparams[2]))*(4**len(fileparams[2]))
            # set number of shots
            newcircuit.attrs['shots'] = shots
            # add raw qasm code to group
            newcircuit.create_dataset('raw_qasm', data=np.string_(qasm))

            # define preparation basis for the set {|0>, |1>, |+>, |i+>}
            basisset = np.asarray(
                ['', 'x q[{0}];', 'h q[{0}];', 'h q[{0}];\ns q[{0}];'])

            # generate all circuit preparation stages and for each all
            # tomography permutations
            for prepcircuit, prep in TCircGen(fileparams[2], basisset):
                for statecircuit, state in TCircGen(fileparams[2]):
                    # create data groups with appropriate formatting
                    grpname = filename[
                        :-5] + '_{}:{}:{}:{}:{}:{}_'.format(*gmtime()[0:6]) + prep + '_' + state
                    tomogcirc = newcircuit.create_group(grpname)
                    # add status attribute: 0 for ready, 1 for uploaded, 2 for
                    # complete
                    tomogcirc.attrs['status'] = int(0)
                    # add idExecution code (null for now)
                    tomogcirc.attrs['executionId'] = np.string_('Null')
                    # insert preparation string
                    prepqasm = qasm.replace(
                        '//PROCESS', prepcircuit + '\n')
                    # add specific qasm code
                    tomogcirc.create_dataset(
                        'qasm', data=np.string_(prepqasm + statecircuit))

        # run circuit as is
        elif fileparams[1] == 'none':
            newcircuit = archive.create_group(dataname)

            # add appropriate attributes

            # device type to run circuits on
            newcircuit.attrs['device'] = np.string_(fileparams[0])
            # type of tomography that is being performed (state or process)
            newcircuit.attrs['tomography'] = np.string_(fileparams[1])
            # target qubits to measure
            newcircuit.attrs['qubits'] = fileparams[2]
            # number of completed circuits (results retrieved)
            newcircuit.attrs['complete'] = 0
            # total number of circuits to be run
            newcircuit.attrs['total'] = 1
            # set number of shots
            newcircuit.attrs['shots'] = shots
            # add raw qasm code to group
            newcircuit.create_dataset('raw_qasm', data=np.string_(qasm))

            # create sub group for circuit
            grpname = filename[:-5] + \
                '_{}:{}:{}:{}:{}:{}_'.format(*gmtime()[0:6])
            singlecirc = newcircuit.create_group(grpname)
            # add status attribute: 0 for ready, 1 for uploaded, 2 for
            # complete
            singlecirc.attrs['status'] = int(0)
            # add idExecution code (null for now)
            singlecirc.attrs['executionId'] = np.string_('Null')
            # create group specific qasm dataset
            singlecirc.create_dataset('qasm', data=np.string_(qasm))

            return
        else:
            archive.flush()
            archive.close()
            print(
                '------- UNRECOGNISED TOMOGRPAHY PARAMETER: \'{}\' ABORTING -------'.format(fileparams[1]))
            exit()

    except ValueError as e:
        # catch group construction error and close archive file with grace
        print(e)
        archive.flush()
        archive.close()
        exit()


# adds target qasm files to archive in batch files from txt file
def batchparse(archive, shots, queuefile='filequeue.txt'):
    # path to experiment data file
    if os.path.isdir('C:\\Users\\Helios'):
        queuepath = 'C:\\Users\\Helios\\Documents\\Projects\\Uni\\2017\\markovpy'
    else:
        queuepath = 'C:\\Users\\Joshua\\Documents\\Projects\\Uni\\2017\\Research\\markovpy'

    # check that target file exists
    path = queuepath + '\\' + queuefile
    if os.path.isfile(path):
        # open file and parse all awaiting circuits
        with open(path, 'r') as queue:
            # extract lines and drop blank lines
            lines = (line.rstrip() for line in queue)
            lines = list(line for line in lines if line)

            # lines that are invalid are rewritten to list queue for user
            rewrite = []
            # parse each new circuit file with parameters
            for index, item in enumerate(lines):
                circuit = item.split(';')
                # check for correct input formatting
                if len(circuit) == 4:
                    file = circuit[0]
                    fileparams = [circuit[1], circuit[
                        2], list(map(int, circuit[3]))]
                    archiveadd(archive, file, fileparams, shots)
                else:
                    rewrite.append(index)
                    print(
                        '------- CIRCUIT {} IN QUEUE FILE HAS INCORRECT NUMBER OF ARGUMENTS: IGNORING -------'.format(index))

        # write all lines that were not parsed back to file
        lines = [lines[i] for i in rewrite]
        with open(path, 'w') as queue:
            for line in lines:
                queue.write(line + '\n')

    else:
        print('------- CIRCUIT BATCH FILE \'{}\' NOT FOUND IN WORKING DIRECTORY: RECOMMEND CREATION -------'.format(path))


# handles experiment interfacing and data arbitration


def experimentrun(archive, interface, filename=None, filepath=None, tomography="ML", shots=1024):
    if qInterface.backend_status('ibmqx2')['available']:
        print('------- IBM QUANTUM COMPUTER IS ONLINE -------')
        ibmOnline = True
    else:
        print('------- IBM QUANTUM COMPUTER IS OFFLINE: UNABLE TO RUN ON EXPERIMENTAL HARDWARE -------')
        ibmOnline = False

    # add new circuits to archive file if any exist
    batchparse(archive, shots)

    # run experiments remotely
    jobcodes = []
    for circuit in archive:
        ccirc = archive[circuit]

        # check if experiment has any circuits awaiting completion for current
        # circuit
        if ccirc.attrs['complete'] < ccirc.attrs['total']:
            device = ccirc.attrs['device'].decode('ascii')

            # iterate over archive circuits, handling status as needed:
            # status = 0: upload file
            # status = 1: retrieve data
            # status = 2: do nothing
            for perm in ccirc:
                # skip raw_qasm, could put in attributes but.....size
                # guidelines
                if perm == 'raw_qasm':
                    continue
                # add files with status 0 (ready) to job queue if qc online or
                # running on simulator
                elif ccirc[perm].attrs['status'] == 0 and (ibmOnline or device == 'simulator'):
                    jobcodes.append(
                        {'qasm': ccirc[perm]['qasm'][()].decode('ascii')})
                # check if past jobs are finished
                elif ccirc[perm].attrs['status'] == 1:
                    # get execution result from IBM
                    result = qInterface.get_result_from_execution(
                        ccirc[perm].attrs['executionId'].decode('ascii'))
                    # check for null result (computation not done yet)
                    if any(result):
                        # check if previous result exists and if so, average measurement result
                        if 'values' in ccirc[perm]:

                            if len(np.asarray(ccirc[perm]['values'][()])) != len(np.asarray(result['measure']['values'])):
                                print('------- RETURN VALUE DIMENSION MISMATCH for {}: ABORTING -------'.format(perm))
                                archive.flush()
                                archive.close()
                                exit()
                            else:
                                newval = 0.5*np.asarray(ccirc[perm]['values'][()]) + 0.5*np.asarray(result['measure']['values'])
                                ccirc[perm].__delitem__('values')
                                ccirc[perm].create_dataset('values', data=newval)

                        else:
                            # add label data to archive file
                            ccirc[perm].create_dataset('labels', data=np.array(
                                result['measure']['labels'], dtype=object), dtype=h5py.special_dtype(vlen=str))
                            # add measured probabilities to archive file
                            ccirc[perm].create_dataset('values', data=result['measure']['values'])

                        # set circuit permutation status to complete
                        ccirc[perm].attrs['status'] = 2
                        # increment number of completed circuits in parent
                        # group
                        ccirc.attrs['complete'] += 1

                # upload jobs in batches of 50
                if len(jobcodes) == 50:
                    if interface.get_my_credits()['remaining'] >= 5 or device == 'simulator':
                        uploads = interface.run_job(
                            jobcodes, backend=device, shots=shots, max_credits=100)
                        codedist(archive, uploads['qasms'])
                    else:
                        print(
                            '------- OUT OF CREDITS TO RUN EXPERIMENT: AWAITING RENEWAL -------')
                    jobcodes = []

    # empty job list if any circuits remain
    if len(jobcodes) > 0:
        if interface.get_my_credits()['remaining'] >= 5 or device == 'simulator':
            uploads = interface.run_job(
                jobcodes, backend=device, shots=shots, max_credits=100)
            codedist(archive, uploads['qasms'])
        else:
            print('------- OUT OF CREDITS TO RUN EXPERIMENT: AWAITING RENEWAL -------')
        jobcodes = []

    # perform state tomography using defined method on all completed circuits
    # instantiate matlab engine
    mateng = matlab.engine.start_matlab()
    for circuit in archive:
        # check if circuit is complete and is compliant with state
        # tomography
        if (archive[circuit].attrs['complete'] == archive[circuit].attrs['total']):
            path = '/' + circuit + '/tomography_' + tomography
            # check if tomography already performed
            if path not in archive:
                    # perform state tomography
                if archive[circuit].attrs['tomography'].decode('ascii') == 'state':
                    print(
                        '------- PERFORMING STATE TOMOGRAPHY ON CIRCUIT: {} -------'.format(circuit))
                    # perform state reconstruction using appropriate method
                    density = getattr(itm, 'statetomography_' +
                                      tomography)(archive, circuit, mateng, goal='state')
                    # add reconstructed operator to archive
                    archive[circuit].create_dataset(
                        'tomography_' + tomography, data=density)
                    # flush archive to be safe
                    archive.flush()

                # perform state tomography on all preperation states, storing
                # reconstructed states in new group
                elif archive[circuit].attrs['tomography'].decode('ascii') == 'process' and tomography == 'ML':
                    # perform state reconstruction using appropriate method
                    state_set = getattr(itm, 'statetomography_ML')(
                        archive, circuit, mateng, goal='process')
                    newcircuit = archive[circuit].create_group('tomography_ML')
                    for pair in state_set:
                        # add reconstructed operator to archive
                        archive[circuit]['tomography_ML'].create_dataset(
                            'tomography_ML_' + pair[0], data=np.asarray(pair[1]))
                        # flush archive to be safe
                        archive.flush()

            # perform process tomography on all completed and reconstructed circuits
            # check if process tomography has already been performed
            path = '/' + circuit + '/Data_Group'
            if path not in archive:
                if archive[circuit].attrs['tomography'].decode('ascii') == 'process' and archive[circuit].get('tomography_ML', getclass=True) == h5py._hl.group.Group:
                    print(
                        '------- PERFORMING PROCESS TOMOGRAPHY ON CIRCUIT: {} -------'.format(circuit))
                    # iterate over computed density matrices
                    chi, opbasis, prepbasis = itm.processtomography(
                        archive, circuit, mateng)
                    # compute kraus operators
                    kraus = itm.mapcompute(chi, opbasis)
                    # compute Choi matrix representation
                    choi = itm.kraus2choi(kraus)
                    # create data group
                    datagroup = archive[circuit].create_group('Data_Group')
                    # add computed chi matrix to archive
                    datagroup.create_dataset('Process_matrix', data=chi)
                    # add operator sum representation
                    datagroup.create_dataset('Kraus_set', data=kraus)
                    # add Choi matrix
                    datagroup.create_dataset('Choi_matrix', data=choi)
                    # add operator basis used to compute above
                    datagroup.create_dataset(
                        'Operator_basis_set', data=opbasis)
                    # add preparation basis used to compute above
                    datagroup.create_dataset(
                        'Preperation_basis_set', data=prepbasis)
                    # flush archive to be safe
                    archive.flush()


#--------------------------------------------------------------------------
# TO BE IMPLEMENTED
#--------------------------------------------------------------------------


#--------------------------------------------------------------------------
# CURRENTLY WORKING ON
#--------------------------------------------------------------------------


if __name__ == '__main__':
    # construct API object for interfacing with the IBM QC
    qInterface = IBMQuantumExperience(token, config)
    # hdf5 data storage file
    with h5py.File(archivepath, 'a') as archive:
        # archive.__delitem__('hadamardq0_simulator')
        try:
            experimentrun(archive, qInterface, shots=8192)
            #itm.densityplot(archive, '5QSTentangle_ibmqx2')
        except Exception as e:
            print(e)
            archive.flush()
            archive.close()
            exit()
        archive.flush()
        archive.close()

    #from gitmanage import gitcommit
    #gitcommit(os.getcwd())
