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
from EZHist import initHist, pushToHist 

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
nHits = 0

def runmaster(args,nClients):

    hd = hitdata()

    inith5(args)
    initHist('DropSize',100,0,2000)
    initHist('TOFhits',1000,-2,2)

    while nClients > 0:
        # Remove client if the run ended
        if hd.recv():
#            print "end run called"
            nClients -= 1
        else:
#            print "Calling plot?"
            plot(hd)
            comp = hd.myobj['comp']
            size=(comp['drop']['a']+comp['drop']['b'])/2.0
            pushToHist('DropSize',size)
            pushToHist('TOFhits',comp['tofsum'])
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

    comp=hd.myobj['comp']
  
    fids[len(fids)-1]=comp['et'].fiducial()
    
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
    global h5out
    fname=args.exprun+'.h5'
    h5out = h5py.File(fname, 'w')

def writeh5(hd):
    global h5out, nHits

    comp = hd.myobj['comp']
    nHits += 1
    hitN = str(nHits)
    h5out[hitN + '/seconds'] = comp['et'].seconds()
    h5out[hitN + '/nanoseconds'] = comp['et'].nanoseconds()
    h5out[hitN + '/fiducial'] = comp['et'].fiducial()    

    h5out[hitN + '/image'] = hd.myorig

    h5out[hitN + '/TOF'] = comp['tof']
    h5out[hitN + '/TOFAxis'] = comp['tofAxis']

    h5out[hitN + '/Dropfit/fitImage'] = hd.myfit
    h5out[hitN + '/Dropfit/a'] = comp['drop']['a']
    h5out[hitN + '/Dropfit/b'] = comp['drop']['b']
    h5out[hitN + '/Dropfit/phi'] = comp['drop']['phi']
    h5out[hitN + '/Dropfit/theta'] = comp['drop']['theta']
    h5out[hitN + '/Dropfit/peakPos'] = comp['drop']['peakPos']
    h5out[hitN + '/Dropfit/peakHeights'] = comp['drop']['peakH']
    h5out[hitN + '/Dropfit/fitFunc'] = comp['drop']['ovalFunc']
    h5out[hitN + '/Dropfit/reducedResidual'] = comp['drop']['reducedRes']

    for name in comp['epics']:
        h5out[hitN + '/epics/' + name] = comp['epics'][name]

def closeh5():
    global h5out
    summaryGroup = h5out.create_group('Summary')
    summaryGroup.create_dataset('message', data = "End_of_run")
    h5out.close()
