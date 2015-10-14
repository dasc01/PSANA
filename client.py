
#run this with: mpirun -n 2 python mpiGather.py
from psana import *
import numpy as np
from hitdata import hitdata 

from fitDroplet import dofitting
from tofHitFind import tofhitfind 

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

def runclient(args):
    hd=hitdata()
    hf=tofhitfind() 

    ds = DataSource(args.exprun+':smd')
    det1 = Detector('pnccdFront',ds.env())
    det2 = Detector('ACQ4',ds.env())
    epics = ds.env().epicsStore()
    thresh = 0.40
    
#    print "rank, size = ", rank, size
    
    for nevent,evt in enumerate(ds.events()):
        if nevent%(size-1)==rank-1:  # different ranks look at different events
#           print rank," processing event", nevent
            eid = evt.get(EventId)
            sec = eid.time()[0]
            nsec = eid.time()[1]
            fid = eid.fiducials()
            et = EventTime(int((sec<<32)|nsec),fid)

            tof, tofAxis = det2.raw(evt)
            
	    if (hf.hitfind(tof[0],tofAxis[0],thresh)):
                epdict=dict()
                for name in epics.names():
                    epdict[name]=epics.value(name)

                img = det1.image(evt)
                obj = dofitting(img)  #call fitting routine
                comp = {'et':et , 'tof':tof[0], 'tofAxis':tofAxis[0], 'epics':epdict, 'drop':obj['drop']}
                hd.send(comp, obj['orig'], obj['fit'])	 
        if nevent == args.noe : break

    hd.endrun()	
    
