#!/usr/bin/python

import ccdproc

modeling = True
if modeling:
    from astropy.modeling import models
from astropy.io import fits

# This works on one multi ext image, not a list, and uses biassec and trimsec
def oscan_trim_file(fname):
  hdulist = fits.open(fname)
  nhdus = len(hdulist)
  if nhdus > 1:
      istart = 1
  else:
      istart = 0
  # loop from first-data to last HDU.
  for i in range(nhdus)[istart:] :
     hdulist = fits.open(fname)
     data1 = ccdproc.CCDData(hdulist[i].data, unit="adu")
     data1.header = hdulist[i].header
     # What happens if file is already overscan-subtracted?
     if modeling:
        oscan1 = ccdproc.subtract_overscan(data1, fits_section=data1.header['BIASSEC'], add_keyword={'overscan': True, 'calstat': 'O'}, model=models.Polynomial1D(1))
     else:
        oscan1 = ccdproc.subtract_overscan(data1, fits_section=data1.header['BIASSEC'], add_keyword={'overscan': True, 'calstat': 'O'}, model=None)
        
     trim1 = ccdproc.trim_image(oscan1, fits_section=oscan1.header['TRIMSEC'], add_keyword={'trimmed': True, 'calstat': 'OT'})
     fits.update(fname, trim1.data, header=trim1.header, ext=i)
  hdulist.close()
  return





