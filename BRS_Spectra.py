import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from photutils.segmentation import SegmentationImage
h = fits.open('DATACUBE_UDF-MOSAIC.fits')
H = h[1].data
h2 = fits.open('MUSE_Segmentation.fits')
h2.info
H2 = h2[0].data

# Find pixels with assigned values in segm
def pixel_mask(Value):
    shape = (945, 947)
    Empty = np.empty(shape)
    for i in range(0, len(H2)):
        Row = H2[i]
        for i2 in range(0, len(Row)):
            Pixel = Row[i2]     
            if Pixel == Value:
                Empty[i2][i]=True
            else:
                Empty[i2][i]=False
    return Empty                    

# Find/extract said pixels in datacube
def pixel_finder(Value):
    Mask = pixel_mask(Value)
    Intensity = []
    for i in range (0, len(H)):
        Slice = H[i]
        int = 0
        for i2 in range(0,len(Mask)):
            Row = Mask[i2]
            Row2 = Slice[i2]
            for i3 in range(0, len(Row)):
                Pixel = Row[i3]
                Pixel2 = Row2[i3]
                if Pixel == True:
                    int = int + Pixel2
                else:
                    int = int
        Intensity = Intensity + [int]
    return Intensity

Wavelength = np.linspace(4750, 9350, 3681)
Intensity_Blue = pixel_finder(774)
Intensity_Red = pixel_finder(729)
Intensity_Star = pixel_finder(1041)

Blue = plt.plot(Wavelength, Intensity_Blue)
Red = plt.plot(Wavelength, Intensity_Red)
Star = plt.plot(Wavelength, Intensity_Star)

Blue.savefig("Blue_Galaxy_Spectrum.pdf")
Red.savefig("Red_Galaxy_Spectrum.pdf")
Star.savefig("Star_Spectrum.pdf")


# Keep each value set separate
# Sum all pixels in a value set together for each slice
# Place combined pixel values across all slices into spectrum graph