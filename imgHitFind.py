import numpy as np


class imghitfind(object):

    sum=0.0
    count=0
    maty=[]
    nav=0

    def __init__(self):
        pass

    def ihitfind(self,img,thresh):
        if(self.count == 0):
            self.maty=np.zeros(img.shape)

        if(self.count < 60):
            self.maty = self.maty + img
            self.nav = self.nav + 1

	self.count += 1 

	self.sum = np.sum(img - self.maty/self.nav)
        print "img:", self.sum, thresh
#        return 1
        return (self.sum > thresh)
