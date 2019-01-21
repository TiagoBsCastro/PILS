from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from glob import glob

for i in range(4):

	kappa = fits.getdata("Alice_PLC/mr_"+str(i).zfill(2)+"/1.0/kappa.fits")
	gamma = fits.getdata("Alice_PLC/mr_"+str(i).zfill(2)+"/1.0/gamma.fits")
	mu = 1.0/((1-kappa)**2-gamma**2)    
	print kappa.mean(), mu.mean()
	print np.average(kappa.flatten(), weights=1.0/mu.flatten()), np.average(kappa.flatten()**2, weights=1.0/mu.flatten())

plt.imshow(kappa)
plt.show()
