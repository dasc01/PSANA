
#run this with: mpirun -n 2 python mpiGather.py
from psana import *
import numpy as np
from hitdata import hitdata 

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

def runclient(args):
    hd=hitdata()
     
    ds = DataSource(args.exprun+':smd')
    det1 = Detector('pnccdFront',ds.env())
    det2 = Detector('ACQ4',ds.env())
    
#    print "rank, size = ", rank, size
    
    for nevent,evt in enumerate(ds.events()):
        if nevent%(size-1)==rank-1:  # different ranks look at different events
#            print rank," processing event", nevent
            img = det1.image(evt)
            intensity = img.sum()
	    if ((nevent)%2 == 0):
	       hd.send(intensity,img)	 
        if nevent == args.noe : break

    hd.endrun()	
    
