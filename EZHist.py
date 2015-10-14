from pypsalg.Histogram import hist1d
from psmon import publish
from psmon.plots import Hist
import numpy as np

histList = dict()

def initHist(histName, nBin, xLow, xHigh):
    global histList
    histList[histName] = hist1d(nBin, xLow, xHigh)

def pushToHist(histName, n):
    global histList
    histList[histName].fill(n)
    xax = np.arange(histList[histName].xaxis.low, histList[histName].xaxis.high+histList[histName].xaxis.binsize, histList[histName].xaxis.binsize)
    histplot = Hist(0, histName, xax, histList[histName].data)
    publish.send(histName, histplot)
