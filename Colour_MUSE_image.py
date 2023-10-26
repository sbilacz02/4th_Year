import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
h = fits.open('DATACUBE_UDF-MOSAIC.fits')
h.info()
H = h[1].data

R = H[0:1226]
G = H[1227:2453]
B = H[2454:3680]

def RGB_file(Colour):
    shape = (947,945)
    Empty = np.empty(shape)
    for i in range(0, len(Colour)):
        Empty = Empty + np.array(Colour[i])

R_image = RGB_file(R)
G_image = RGB_file(G)
B_image = RGB_file(B)

hdu_R = fits.PrimaryHDU(R_image)
hdul_R = fits.HDUList([hdu_R])
hdul_R.writeto('R_Image.fits')

hdu_G = fits.PrimaryHDU(G_image)
hdul_G = fits.HDUList([hdu_G])
hdul_G.writeto('G_Image.fits') 

hdu_B = fits.PrimaryHDU(B_image)
hdul_B = fits.HDUList([hdu_B])
hdul_B.writeto('B_Image.fits') 