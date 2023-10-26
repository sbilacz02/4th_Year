import numpy as np
from astropy.io import fits
h = fits.open('DATACUBE_UDF-MOSAIC.fits')
h.info()
H = h[1].data
