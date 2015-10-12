from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

from psmon import publish
from psmon.plots import Image,MultiPlot
import h5py
import numpy as np
from hitdata import hitdata 

# might need classes/routines to:
# - plot
# - remember the history
# - save the history to h5
npanel=5
figs=[None]*npanel

def runmaster(nClients):
    hd = hitdata()
    
    while nClients > 0:
        # Remove client if the run ended
        if hd.recv():
#            print "end run called"
            nClients -= 1
        else:
#            print "Calling plot?"
            plot(hd)

def plot(hd):
    print hd.myimg.shape
    for j in range(0,len(figs)-1):
	figs[j]=figs[j+1]
    figs[len(figs)-1]=hd.myimg   
    multi = MultiPlot(1, 'Some Plots')
    for j in range(0,len(figs)):
        if figs[j] is not None :
      	    plotimg = Image(j,"CsPad",figs[j])
            multi.add(plotimg)
    publish.send('MULTI', multi)
#    time.sleep(0.4)
