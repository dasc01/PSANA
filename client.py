#run this with: mpirun -n 2 python mpiGather.py
from psana import *
import numpy as np
from hitdata import hitdata 

from fitDroplet import dofitting
from tofHitFind import tofhitfind
from imgHitFind import imghitfind

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

def runclient(args, mask):
    hd=hitdata()
    hf=tofhitfind() 
    ihf=imghitfind()
    
#    ds = DataSource(args.exprun+':smd')
    ds = DataSource(args.exprun+':smd'+':dir=/reg/d/ffb/amo/amoj5415/xtc:live')
    det1 = Detector('pnccdFront',ds.env())
    det2 = Detector('ACQ4',ds.env())
    epics = ds.env().epicsStore()
    thresh = 0.0
    ithr   = 1E8
    
#    print "rank, size = ", rank, size
    
    for nevent,evt in enumerate(ds.events()):
        if nevent%(size-1)==rank-1:  # different ranks look at different events
            try:
                print rank," processing event", nevent
                eid = evt.get(EventId)
                sec = eid.time()[0]
                nsec = eid.time()[1]
                fid = eid.fiducials()
                et = EventTime(int((sec<<32)|nsec),fid)

		tof = det2.raw(evt)
		if tof is not None: 
		    tof, tofAxis = det2.raw(evt)
		else:
		    print nevent, " is faulty"
		    continue
            
                if (hf.hitfind(tof[0],tofAxis[0],thresh)):
                    epdict=dict()
                    for name in epics.names():
                        epdict[name]=epics.value(name)

                    img = det1.image(evt)
		    if(np.random.random() > 0.99):   #ihf.ihitfind(img,ithr)):
                        obj = dofitting(img, mask)  #call fitting routine
                        comp = { 'et':et , 'tof':tof[0], 'tofAxis':tofAxis[0], 'epics':epdict, 'drop':obj['drop'], 'tofsum':hf.sum , 'imgsum':ihf.sum }
                        hd.send(comp, obj['orig'], obj['fit'])	 

                if nevent == args.noe : break

            except:

                print "Exception on rank: ", rank, " event ", nevent
                continue


    hd.endrun()	
    
