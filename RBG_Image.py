import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils.segmentation import SegmentationImage
h = fits.open('DATACUBE_UDF-MOSAIC.fits')
h.info()
H = h[1].data

Super_Empty = []
for i in range(0, len(H2)):
    Section = H2[i]
    Empty = []
    for i2 in range(0, len(Section)):
        Empty = Empty + int(Section[i2])
    Super_Empty = Super_Empty + Empty

R = H[2454:3680]
G = H[1227:2453]
B = H[0:1226]

def RGB_file(Colour):
    shape = (947,945)
    Empty = np.empty(shape)
    for i in range(0, 1226):
        Empty = Empty + np.array(Colour[i])
    return Empty

R_image = RGB_file(R)
G_image = RGB_file(G)
B_image = RGB_file(B)

def getstats(img,ff):
    segm = SegmentationImage(img)
    gd = np.where((img != 0) & (img != np.nan))
    print('there are ',len(gd[0]),' elements')
    arr = img[gd]
    arr = sorted(arr)
    n = len(arr)
    print('array is ',n,' elements')
    i = round(ff*n)
    vmax = arr[i]
    print(ff,' signal range value is ',vmax)
    
    print('making mask')
    mask = segm.make_source_mask(img, nsigma=2, npixels=5, dilate_size=11)
    print('calculating stats')
    vmean, vmedian, vstd = sigma_clipped_stats(img, sigma=3.0, mask=mask,mask_value=0.)
    print('mean: ',vmean)
    print('median: ',vmedian)
    print('sigma: ',vstd)
    return vmean,vmedian,vstd,vmax

def mkcol(b,v,r,ff,gamma,xlo,xhi,ylo,yhi):

#    xlo = 300#4000#1700
#    xhi = 1500#6000#3100
#    ylo = 300#15000#6000
#    yhi = 1500#17000#7500

    tmpb = b[ylo:yhi,xlo:xhi]
    tmpv = v[ylo:yhi,xlo:xhi]
    tmpr = r[ylo:yhi,xlo:xhi]

    #tmpb = b#[ylo:yhi,xlo:xhi]
    #tmpv = v#[ylo:yhi,xlo:xhi]
    #tmpr = r#[ylo:yhi,xlo:xhi]

    #tmpb = nan2d(tmpb,0)
    #tmpv = nan2d(tmpv,0)
    #tmpr = nan2d(tmpr,0)

    #print(np.sum(tmpb))
    #print(np.sum(tmpv))
    #print(np.sum(tmpr))
    
    bmean,bmedian,bstd,bmax = getstats(tmpb,ff)
    vmean,vmedian,vstd,vmax = getstats(tmpv,ff)
    rmean,rmedian,rstd,rmax = getstats(tmpr,ff)

    print('rescaling...')
    bmin = bmean
    vmin = vmean
    rmin = rmean

    gdb = np.where((b != 0) & (b!=np.nan))
    gdv = np.where((v != 0) & (v!=np.nan))
    gdr = np.where((r != 0) & (r!=np.nan))

    b[gdb] = (b[gdb]-bmin)/(bmax-bmin)
    v[gdv] = (v[gdv]-vmin)/(vmax-vmin)
    r[gdr] = (r[gdr]-rmin)/(rmax-rmin)

    lo = 0.
    hi = 1.

    bad = np.where(b <= lo)
    b[bad]=0.
    bad = np.where(b >= hi)
    b[bad]=1.

    bad = np.where(v <= lo)
    v[bad]=0
    bad = np.where(v >= hi)
    v[bad]=1.

    bad = np.where(r <= lo)
    r[bad]=0
    bad = np.where(r >= hi)
    r[bad]=1

    b = b**gamma
    v = v**gamma
    r = r**gamma

    print('writing to array')
#    b = b*254.
#    v = v*254.
#    r = r*254.

    sz = b.shape
    #print(sz[1],sz[0])
    
    col = np.zeros((sz[0],sz[1],3))
    col[:,:,0] = r
    col[:,:,1] = v
    col[:,:,2] = b
    print('mkcol complete.')
    return col

RGB = mkcol(B_image, G_image, R_image, 0.98, 0.8)

hdu = fits.PrimaryHDU(RGB)
hdul = fits.HDUList([hdu])
hdul.writeto('RGB_Image.fits', overwrite=True)
