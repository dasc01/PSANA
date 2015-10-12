from psana import *
import numpy as np
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()



class hitdata(object):

    myobj={}
    myimg=[]

    def __init__(self):
        pass

    def endrun(self):
        obj={'endrun':True}
#        print rank,"is ending run"
        comm.send(obj,dest=0,tag=rank)
#        print rank,"has ended run"
        

    def send(self, et, img):
#        print rank,"is about to send"
        obj={'et':et,'shape':img.shape,'endrun':False}
        comm.send(obj,dest=0,tag=rank)
        comm.Send([img,MPI.DOUBLE],dest=0,tag=rank)
#        print rank, "has sent"
 

    def recv(self):
        status=MPI.Status()       
#        print rank,"is waiting to recv"
        self.myobj=comm.recv(source=MPI.ANY_SOURCE,tag=MPI.ANY_TAG,status=status)
        recvRank = status.Get_source()
        if self.myobj['endrun'] == False:
           self.myimg=np.empty(self.myobj['shape'],dtype=np.float32)
           comm.Recv(self.myimg,source=recvRank,tag=MPI.ANY_TAG)
#        print rank,"recvd from", recvRank, self.myobj['endrun']
        return (self.myobj['endrun'])



