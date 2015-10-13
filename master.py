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
ftop=[None]*npanel
fbot=[None]*npanel
fmid=[None]*npanel
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

    for j in range(0,len(ftop)-1):
	ftop[j]=ftop[j+1]
	fmid[j]=fmid[j+1]
        fids[j]=fids[j+1]

    ftop[len(ftop)-1]=hd.myorig
    fmid[len(fmid)-1]=hd.myfit  
    fids[len(fids)-1]=hd.myobj['et'].fiducial()
    
    multop = MultiPlot(1, 'Original Image')
    mulmid = MultiPlot(1, 'Fitted   Image')

    for j in range(0,len(ftop)):
        if ftop[j] is not None :
      	    plottop = Image(fids[j],"Original",ftop[j])
      	    plotmid = Image(fids[j],"Fitted",fmid[j])
            multop.add(plottop)
            mulmid.add(plotmid)

    publish.send('ORIG', multop)
    publish.send('FIT', mulmid)
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
