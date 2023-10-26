import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
h = fits.open('DATACUBE_UDF-MOSAIC.fits')
H = h[1].data

Wavelength = np.linspace(wI, wF, 3680)
Intensity = []
for i in range(0, len(H)):
    Target = H[i]
    Intensity = Intensity + Target[x][y]

plt.plot(Wavelength, Intensity)