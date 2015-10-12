from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

from psmon import publish
from psmon.plots import Image,MultiPlot
import h5py
import numpy as np
from hitdata import hitdata 
import time

# might need classes/routines to:
# - plot
# - remember the history
# - save the history to h5

npanel=5
figs=[None]*npanel
fids=[None]*npanel

#Set up hdf5 stuff
h5out = None
eventDataGroup=None
evtSecDs=None
evtNanoDs=None
nextDsIdx=-1

def runmaster(args,nClients):

    hd = hitdata()

    inith5(args)

    while nClients > 0:
        # Remove client if the run ended
        if hd.recv():
#            print "end run called"
            nClients -= 1
        else:
#            print "Calling plot?"
            plot(hd)
            writeh5(hd)

    closeh5()

#PLOT function-------------------------
def plot(hd):

    for j in range(0,len(figs)-1):
	figs[j]=figs[j+1]
        fids[j]=fids[j+1]

    figs[len(figs)-1]=hd.myimg
    fids[len(fids)-1]=hd.myobj['et'].fiducial()
    
    multi = MultiPlot(1, 'Some Plots')

    for j in range(0,len(figs)):
        if figs[j] is not None :
      	    plotimg = Image(fids[j],"CsPad",figs[j])
            multi.add(plotimg)

    publish.send('MULTI', multi)
#    time.sleep(1.0)

#HDF5 functions--------------------------
def inith5(args):
    global h5out, eventDataGroup, evtSecDs, evtNanoDs
    fname=args.exprun+'.h5'
    h5out = h5py.File(fname, 'w')
    eventDataGroup = h5out.create_group('EventData')
    evtSecDs = eventDataGroup.create_dataset('event_seconds',(0,), dtype='i4', chunks=True, maxshape=(None,))
    evtNanoDs = eventDataGroup.create_dataset('event_nanoseconds',(0,), dtype='i4', chunks=True, maxshape=(None,))


def writeh5(hd):
    global evtSecDs, evtNanoDs, nextDsIdx

    nextDsIdx += 1
 
    evtSecDs.resize((nextDsIdx+1,))
    evtNanoDs.resize((nextDsIdx+1,))
 
    evtSecDs[nextDsIdx] = hd.myobj['et'].seconds()
    evtNanoDs[nextDsIdx] = hd.myobj['et'].nanoseconds()

def closeh5():
    global h5out
    summaryGroup = h5out.create_group('Summary')
    summaryGroup.create_dataset('message', data = "End_of_run")
    h5out.close()
