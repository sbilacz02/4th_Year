import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
h = fits.open('DATACUBE_UDF-MOSAIC.fits')
H = h[1].data

image = H[458]

plt.imshow(image, cmap = 'gray')

hdu = fits.PrimaryHDU(image)
hdul = fits.HDUList([hdu])
hdul.writeto('MUSE_Image.fits')