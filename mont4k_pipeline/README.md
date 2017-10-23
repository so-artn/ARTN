# Mont4K pipeline

This is a python / astropy / ccdproc pipeline for data from the 
Mont4K CCD imager on Steward Observatory's Kuiper 61" telescope.

It currently will take a directory of data, determine which files
are biases / dome flats /sky flats / objects,  do overscan subtraction
and trimming, make a master bias, do bias subtraction, make master
flats, flatfield the objects with dome flats, and merge the 2 FITS
extensions from the 2 readouts of the CCD into a single FITS extension.
It uses the FILTER keyword to associate objects with the correct flat.

The only part of this pipeline that is really Mont4K specific is the
merging into a single FITS image and the geometry of the array. Much
of the code could be used to reduce multi-extension FITS images (MEF)
from other instruments. This may be useful for others, since astropy's ccdproc 
doesn't natively support multi-extension FITS files, so the convenience
routines in this pipeline had to be written to run on MEFs.

Possibly to be added:

* sky flat / illumination correction
* fringe correction
* superskies?
* cosmic ray removal / astroscrappy
* astrometry
* photometry of standards?


