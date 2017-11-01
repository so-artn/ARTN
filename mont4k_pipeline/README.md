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

The pipeline requires astropy and ccdproc. If you don't have these
try "pip install ccdproc" at the command line, etc.
It also currently requires the packages msumastro and bottleneck, but these aren't
really being used and could be commented out.

Look at the comments at the beginning for hwo to run it. You can try:

to process all the files in the current directory, try:

```
python
import sys 
sys.path.append('/Users/bjw/stellarview/python')  # Replace the directory with the directory where you put kuiper_reduce.py
import kuiper_reduce as kp
fstruct = kp.getheaders('')
imagelist = []
for image in fstruct:
   imagelist.append(image['filename'])
kp.process_and_merge_list(imagelist)
```

If you read through the code, most of the routines are pretty modular so you can
run individual routines, use only certain stages of the processing, etc.

