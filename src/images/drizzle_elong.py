from astropy.io import fits
import astropy.wcs as wcs
import reproject as rp
import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as ndi
from drizzle import drizzle
from astropy.nddata import CCDData


import sys

root=sys.argv[1]

d1 = fits.open("%s.detector0_000000_stack.fits" % (root))
#d2 = fits.open("test0418._detector0_000002_stack.fits")
d1i_head = d1[0].header.copy()
print(d1[0].header['CDELT1'])
d0 = d1.copy()
rscl=3
d1i_head = d0[0].header
d0[0].data = np.zeros(shape=(int(d1i_head['NAXIS1']*rscl),int(d1i_head['NAXIS2']*rscl)))
d1i_head['CD1_1'] = d1[0].header['CD1_1']/rscl
d1i_head['CD1_2'] = d1[0].header['CD1_2']/rscl
d1i_head['CD2_1'] = d1[0].header['CD2_1']/rscl
d1i_head['CD2_2'] = d1[0].header['CD2_2']/rscl
d1i_head['CDELT1'] = d1[0].header['CDELT1']/rscl
d1i_head['CDELT2'] = d1[0].header['CDELT2']/rscl
d1i_head['CRPIX1'] = d1[0].header['CRPIX1']*rscl
d1i_head['CRPIX2'] = d1[0].header['CRPIX2']*rscl
owcs = wcs.WCS(d1i_head)


print(owcs)

nframes=36

for j in range(2):
    driz1 = drizzle.Drizzle(outwcs=owcs,pixfrac=0.0)

    for i in range(nframes):
        print(i)
        driz1.add_fits_file("%s.detector%d_%06d_stack.fits" % (root,j,i))

    outname = "%s.%d_drizzle.stack.fits" % (root,j)
    driz1.write(outname)
    print(f"write {outname}")
    elong = fits.open(outname)


    driz1 = drizzle.Drizzle(outwcs=owcs,pixfrac=0.0)

    for i in range(nframes):
        print(i)
        driz1.add_fits_file("%s.detector%d_%06d_nestack.fits" % (root,j,i))
                       
    outname = "%s.%d_drizzle.nestack.fits" % (root,j)
    driz1.write(outname)
    print(f"write {outname}")
    nelong = fits.open(outname)

    diff = elong['SCI'].data - nelong['SCI'].data
    fracdiff = diff/nelong['SCI'].data

    difffits = CCDData(diff,unit='adu')
    difffits.write("%s.%d_drizzle.diff.fits" % (root,j),overwrite=True)

    fracdifffits = CCDData(fracdiff,unit='adu')
    fracdifffits.write("%s.%d_drizzle.fracdiff.fits" % (root,j),overwrite=True)

