#!/usr/bin/env python

import pyregion
import scipy.stats
from astropy.io import fits
from astropy import wcs
import numpy as np
import sys


fitsfile=fits.open(sys.argv[1])
d=fitsfile[1].data
for f in sys.argv[2:]:
    fg_ir=pyregion.open(f).as_imagecoord(fitsfile[1].header)
    for n,reg in enumerate(fg_ir):
        fg=pyregion.ShapeList([reg])
        mask=fg.get_mask(hdu=fitsfile[1],shape=np.shape(d))
        pixels=np.sum(mask)
        if (n==0):
            data=np.extract(mask,d)
        else:
            hld=np.extract(mask,d)
            data=np.append(data,hld)
    
    print('%s %s RMS = %g mean = %g npixel = %d' % (f,f.replace('isoregion_','').replace('.reg',''),scipy.stats.nanstd(data),scipy.stats.nanmean(data),len(data)))

fitsfile.close()


