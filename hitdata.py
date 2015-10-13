from psana import *
import numpy as np
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()



class hitdata(object):

    myobj={}
    myorig=[]
    myfit=[]
    mydrop={}
    means=[]

    def __init__(self):
        pass

    def endrun(self):
        obj={'endrun':True}
#        print rank,"is ending run"
        comm.send(obj,dest=0,tag=rank)
#        print rank,"has ended run"
        

    def send(self, et, orig, fit, drop):
#        print rank,"is about to send"
        obj={'et':et,'shape':orig.shape,'endrun':False}
        comm.send(obj,dest=0,tag=rank)
        comm.Send([orig,MPI.DOUBLE],dest=0,tag=rank+1)
        comm.Send([fit,MPI.DOUBLE],dest=0,tag=rank+2)
#        comm.send(drop,dest=0,tag=rank+3)
#        print rank, "has sent"
 

    def recv(self):
        status=MPI.Status()       
#        print rank,"is waiting to recv"
        self.myobj=comm.recv(source=MPI.ANY_SOURCE,tag=MPI.ANY_TAG,status=status)
        recvRank = status.Get_source()
        if self.myobj['endrun'] == False:
           self.myorig=np.empty(self.myobj['shape'],dtype=np.float32)
           comm.Recv(self.myorig,source=recvRank,tag=recvRank+1)
           self.myfit=np.empty(self.myobj['shape'],dtype=np.float32)
           comm.Recv(self.myfit,source=recvRank,tag=recvRank+2)
#           self.mydrop=comm.recv(source=recvRank,tag=recvRank+3)
        return (self.myobj['endrun'])



