#%matplotlib inline

import scipy.special
import scipy.ndimage
import scipy.optimize
import scipy.fftpack
import scipy.signal
import scipy.stats
import numpy as np
import matplotlib.pylab as plt

def simulateImage(a_nm, b_nm, x0, y0, ph0, mask):
    tiny = 1.0e-16
    a = a_nm*1e-6 #turn the nm drops into mm
    b = b_nm*1e-6

    dpix = 0.075 #%pixels in mm
    d = 567; #%centerline distance from droplet to detector = 710 mm
    k = 7.6e6 #%  X-ray wavenumber = 4.25e6 mm-1 = 2pi/lambda
    #k = 7.6e6/3

    imageSize = 512# %N.B we could make this smaller, because we'll mask most out

    #get a coordinate grid
    X, Y = np.meshgrid((np.arange(0,imageSize)-x0)*dpix, (np.arange(0,imageSize)-y0)*dpix) 
    Phi = np.arctan(Y/(X+tiny)) #%angle grid
    R = np.sqrt(X**2+Y**2) #%radius grid

    Reff = np.sqrt(a**2*(np.cos(Phi-ph0))**2+b**2*(np.sin(Phi-ph0))**2) #%r in oval coords
    Rd = np.sqrt(R**2+d**2) #%actual distance from the FEL focus to the spot hit on the detector

    #% Q is crystal momentum
    ThetA = np.arctan(R/d) #%angle from beam axis?
    Q = 2*k*np.sin(ThetA/2) 
    Q[mask<1]=1 #%?????

    Res1=2*np.pi*(k**4*a**2*b**4)
    Res2=(scipy.special.jv(3/2,Reff*Q))**2
    Res3=1/((Reff*Q)*(Reff*Q)*(Reff*Q))
    Res4=1/(Rd*Rd)
    ReS=Res1*Res2*Res3*Res4
    ReS[np.isnan(ReS)]=0
    return mask*ReS#/np.sum(ReS)


def imageToPolar(image, centreX, centreY):
    angularSize = np.size(image,1)/2;
    radialSize = np.size(image,1)/2;

    th = np.linspace(0, 2*np.pi, angularSize)
    r = np.linspace(0, radialSize*2-centreX, radialSize)
    
    R = np.tile(r, [angularSize,1])
    TH = np.tile(th, [radialSize, 1]).T
    [XI, YI] = pol2cart(R, TH)
    
    polarImage = scipy.ndimage.map_coordinates(image, [XI + centreX, YI + centreY], order=1).T  
    polarImage[np.isnan(polarImage)] = 2
    polarImage[polarImage == 0] = 2
    return polarImage

def getSymmetry(dropletImage, mask, x, y):
    #getSymmetry Check the symmetry of a a polar image.
    #  Takes the image, the mask (0||1) and the centre. Returns the difference
    #  of the 1st and 2nd 180 degree segements, for use in a least squares
    #  fit.

    imLength = np.shape(dropletImage)[0]/2
    polarMask = imageToPolar(mask, x, y)
    totalMask = polarMask[:,0:imLength/2] * polarMask[:,imLength/2:]

    polarImage = imageToPolar(dropletImage, x, y)
    return (polarImage[:,0:imLength/2] - polarImage[:,imLength/2:]) * totalMask

def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)

def pol2cart(rho, phi):
    x = rho * np.cos(phi+np.pi)
    y = rho * np.sin(phi+np.pi)
    return(x, y)
    
def estimateSize(peakSpacing):
    #tweak paramters for given run
    wavelength=0.82 #wavelength (nm)
    deltapix=peakSpacing
    detlamm=deltapix*0.075 #pixel size (in mm)
    deltaphi=(detlamm)/565.5 #interaction --> detector (in mm)
    D =wavelength/deltaphi
    return D/2 
        
def findCentre(image, mask):
    #findCentre Gets the centre of a given image. Needs to be within ~5 pixels 
    #to start with 
    #264.6, 275.4 (x,y for test image)
    fun = lambda x: np.ravel(getSymmetry(image, mask, x[0], x[1]))
    res = scipy.optimize.leastsq(fun, [264, 275], xtol=0.001)
    return res[0]

def fitPolarImage(polarDroplet, angularBins):
    tiny=1.0e-16
    #some parameters for pink noise subtraction:
    pinkExponent = 1.0
    pinkCoef = 2.0
    
    theta = np.linspace(0, 2*np.pi, angularBins)    
    
    #rebin
    if (polarDroplet.shape[0] % angularBins):
        raise ValueError('angularBins must be an even divisor of image size')
    polarDropReshaped = np.reshape(polarDroplet, [polarDroplet.shape[0], angularBins, polarDroplet.shape[1]/angularBins])
    polarDropBinned = polarDropReshaped.mean(axis=2)
    polarDropBinned = polarDropBinned - np.mean(polarDropBinned) 
   
    #TAKE FFT
    dropFFT = np.abs(scipy.fftpack.fft(polarDropBinned, axis=0))[0:polarDropBinned.shape[0]/2,:]
    freqAxis = scipy.fftpack.fftfreq(polarDropBinned.shape[0])[0:polarDropBinned.shape[0]/2]
    meanFFT = np.mean(dropFFT, axis=1)
    pinkScale = pinkCoef*(meanFFT[1]/freqAxis.size) * (255**(1-pinkExponent))
    subtractedFFT = dropFFT - np.tile(pinkScale/(freqAxis**pinkExponent + tiny), [angularBins, 1]).T
  
  #  plt.imshow(polarDropBinned)
  #  plt.show()    
  #  plt.imshow(dropFFT)
  #  plt.show()  
  #  plt.plot(freqAxis, subtractedFFT.mean(axis=1))
  #  plt.plot(freqAxis, meanFFT)
    
    #FIND PEAKS IN FFT
    peakPositions = np.zeros(angularBins)
    peakHeights = np.zeros(angularBins)
    for row in range(0,angularBins):
        angleSlice = subtractedFFT[:,row]
        peakPositions[row] = freqAxis[np.argmax(angleSlice)]
        peakHeights[row] = angleSlice[np.argmax(angleSlice)]
    
    #do some outlier removal
    outlierFactor = 2
    peakHeights[peakPositions > scipy.stats.mode(peakPositions)[0]*outlierFactor] = 0
    peakHeights[peakPositions < scipy.stats.mode(peakPositions)[0]/outlierFactor] = 0
    
#    plt.scatter(theta, peakPositions, s=peakHeights*100)
    #plt.ylim([])

    #FIT PEAKS TO SIN
    ovalFunc = lambda x: (x[0]*np.cos(2*theta + x[1])+x[2])
    fitFun = lambda x: (ovalFunc(x)-peakPositions)*(peakHeights)
    res = scipy.optimize.leastsq(fitFun, [0.0, 0.0, np.sum(peakPositions*peakHeights)/np.sum(peakHeights)])#, xtol=0.001)

#    plt.plot(theta, ovalFunc(res[0]), 'r')
#    plt.ylim([res[0][2] - 4*res[0][0], res[0][2] + 4*res[0][0]])
#    plt.show()
    
    residual = np.sum((ovalFunc(res[0]) - peakPositions)**2 / np.sum(peakHeights))
    
    a = estimateSize(1/(res[0][2] + res[0][0]))
    b = estimateSize(1/(res[0][2] - res[0][0]))
    
    if (a<b):
        temp = a
        a = b
        b = temp
        phi = res[0][1]/2#-np.pi/2
    else:
        phi = res[0][1]/2-np.pi/2
    
    print('a = ' + str(a))
    print('b = ' + str(b))
    print('Reduced residual = ' + str(residual))
    return {'a':a, 'b':b ,'phi':phi, 'theta':theta, 'peakPos':peakPositions, 'peakH':peakHeights, 'ovalFunc':ovalFunc(res[0]) , 'reducedRes':residual }


    #APPROXIMATE DROPLET PARAMETERS
    
def cartesianFit(a, b, x0, y0, phi, logDropletImage, mask):
    #DEFINE FIT FUNC
    #DO A FIT!
    fun = lambda x: np.ravel(np.log10(100+5e6*x[3]*simulateImage(x[0], x[1], x0, y0, x[2], mask)) - logDropletImage)
    res = scipy.optimize.leastsq(fun, [a, b, phi, 1], xtol=0.001, full_output=0)
    print(res)
    plt.imshow(np.log10(100+5e6*res[0][3]*simulateImage(res[0][0], res[0][1], x0, y0, res[0][2], mask)))
    plt.colorbar()
    plt.show()    
    plt.imshow(np.reshape(fun(res[0]), [512, 512]))
    plt.colorbar()
    plt.show()

##NB this is just test stuff, in reality we'll pull from XTC
##39 - tricky, #33 low freq\
##/Users/adam/Dropbox/LCLS_OCT15/Code/run_99/logScale#
#means = np.zeros(50)
#for n in range(0,3):
#    scaledLogImage = np.loadtxt('/reg/neh/home/dasc/PSANA/RawData/run_99/logScale' + str(n) +'.txt', delimiter=',')
#    rawMask = np.loadtxt('mask.txt', delimiter=',')
#
#    #we're going to put these images on a 512x512 basis
#
#    newIm = 2*np.ones([512,512]);
#    newMask = np.zeros([512,512]);
#    newIm[56:56+398,56:56+398] = scaledLogImage;
#    newMask[56:56+398,56:56+398] = rawMask;
#    plt.imshow(newIm)
#    plt.colorbar()
#    plt.show()
#    centre = findCentre(newIm, newMask)
#    print('x0 = ' + str(centre[0]))
#    print('y0 = ' + str(centre[1]))
#    polarDroplet = imageToPolar(newIm, centre[0], centre[1])
#    drop = fitPolarImage(polarDroplet, 32)
#    means[n] = drop['meanSize']
#    plt.imshow(np.log10(100+5.0e5*simulateImage(drop['a'], drop['b'], centre[0], centre[1], drop['phi'], newMask)))
#    plt.colorbar()
#    plt.show()
#    #cartesianFit(drop['a'], drop['b'], centre[1], centre[0], drop['phi'], newIm, newMask)
#    
#plt.hist(means)
#plt.show()

def dofitting(img):
     min = np.amin(img)
     scaledLogImage = np.log10(abs(min)+100+img)[0:1024,0:1024]
     newIm = 2*np.zeros([512,512],dtype=np.float32);
     for i in range(512):
         for j in range(512):
             newIm[i,j]=(scaledLogImage[2*i,2*j]+scaledLogImage[2*i+1,2*j]+scaledLogImage[2*i,2*j+1]+scaledLogImage[2*i+1,2*j+1])/4.0         
     rawMask = np.loadtxt('mask.txt', delimiter=',')
     newMask = np.zeros([512,512],dtype=np.float32);
     newMask[56:56+398,56:56+398] = rawMask;
     centre = findCentre(newIm, newMask)
     polarDroplet = imageToPolar(newIm, centre[0], centre[1])
     drop = fitPolarImage(polarDroplet, 32)
     fitImage = simulateImage(drop['a'], drop['b'], centre[0], centre[1], drop['phi'], newMask).astype(np.float32)
 
     return {'orig':img[0:1024,0:1024] , 'fit':fitImage, 'drop':drop}

#     means[n] = drop['meanSize']
#     plt.imshow()
#     plt.colorbar()
#     plt.show()
#     #cartesianFit(drop['a'], drop['b'], centre[1], centre[0], drop['phi'], newIm, newMask)
    
#plt.hist(means)
#plt.show()


