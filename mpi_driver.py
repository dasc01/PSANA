# run these commands in the "amoh1315" directory

# mpirun -n 2 python mpi_driver.py  exp=amoh1315:run=256 -n 200
# in batch:
# bsub -q psanaq -n 2 -o %J.log -a mympi python mpi_driver.py exp=amotut13:run=206

from master import runmaster
from client import runclient

import numpy as np
from psana import *

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
assert size>1, 'At least 2 MPI ranks required'
numClients = size-1

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("exprun", help="psana experiment/run string (e.g. exp=amoj5415:run=43)")
parser.add_argument("-n","--noe",help="number of events, all events=0",default=-1, type=int)

args = parser.parse_args()
mask = np.empty([1024, 1024])

#pull the mask from the 1st event
ds = DataSource(args.exprun+':smd')
det = Detector('pnccdFront',ds.env())
for nevent,evt in enumerate(ds.events()):
    mask = det.image(evt,np.ones((4,1024,1024)))
    break

from psmon import publish
from psmon.plots import Image

mask = np.ones([1024, 1024])

if rank==0:
    maskplot = Image(0,'mask',mask)
    publish.send('mask', maskplot)
    runmaster(args,numClients, mask)
else:
    runclient(args, mask)


MPI.Finalize()
