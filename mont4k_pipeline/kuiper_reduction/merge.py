#!/usr/bin/python


from astropy.io import fits
import ccdproc
import numpy as np

def merge_fitsfd(hdulist):

	data1 = ccdproc.CCDData(hdulist[1].data, unit="adu")
	data1.header = hdulist[1].header
	data2 = ccdproc.CCDData(hdulist[2].data, unit="adu")
	merged = np.concatenate( (data1, np.fliplr(data2) ), axis=1)
	# assume we don't have to change any WCS parameters from ext 1
	hdu_new = fits.PrimaryHDU( merged )
	hdu_new.header = hdulist[0].header
	hdulist_new = fits.HDUList( [hdu_new] )
	# Use basename to keep from breaking if the names in imagelist are
	# full path rather than just a file

	return hdulist_new




