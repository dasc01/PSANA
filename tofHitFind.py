import numpy as np


class tofhitfind(object):

    sum=0.0
    count=0
    maty=[]
    nav=0

    def __init__(self):
        pass

    def hitfind(self,tof,tofAxis,thresh):
        if(self.count == 0):
            self.maty=np.zeros(len(tof))

        if(self.count < 60):
            self.maty = self.maty + tof
            self.nav = self.nav + 1

	self.count += 1 

	self.sum = 0.0
        for i in range(len(tof)):
            self.sum = self.sum + (tof[i] - self.maty[i]/self.nav)
        print "tof:", self.sum, thresh
#        return 1
        return (abs(self.sum) > thresh)
