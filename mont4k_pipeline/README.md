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
really being used and could be commented out.  But for now you should try
"pip install msumastro", "pip install bottleneck", "pip install ccdproc" before running.

Look at the comments at the beginning for how to run it. You can try:

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

The control flow of the overall package in kuiper_reduce.py 
goes through these routines:

- getheaders
- create imagelist
- process_and_merge_list
  - process_list
    - makelists
    - oscan_trim_list(zerolist)
      - oscan_trim_file
        - ccdproc.subtract_overscan
	-  ccdproc.trim_image
    - make_master_bias
    - oscan_trim_list(domelist)
      - oscan_trim_file
        - ccdproc.subtract_overscan
        - ccdproc.trim_image
    - bias_list_by_file(domelist,masterbiasname)
      - bias_file_by_file
        - ccdproc.subtract_bias
    - make_master_flats(domelist, ...)
      - make_master_oneflat
        - combine_list_to_file
          - ccdproc.Combiner
    - oscan_trim_list(skylist)
      - oscan_trim_file
        - ccdproc.subtract_overscan
        - ccdproc.trim_image
    - bias_list_by_file(skylist,masterbiasname)
      - bias_file_by_file
        - ccdproc.subtract_bias
    - make_master_flats(skylist, ...)
      - make_master_oneflat
        - combine_list_to_file
          - ccdproc.Combiner
    - oscan_trim_list(objectlist)
      - oscan_trim_file
        - ccdproc.subtract_overscan
        - ccdproc.trim_image
    - proc_objects_all(objectlist, masterbiasname, flatnames)
      - find_filters_list
      - proc_objects_filter
        - bias_file_by_file(fname, master_bias)
          - ccdproc.subtract_bias
        - flat_file_by_file(fname, master_flat)
          - ccdproc.flat_correct
          - renormalize_by_flat *This should be either used or deprecated since ccdproc.flat_correct can now take a flat scale as input*

  - merge_m4k_list (specific to the Mont4K 2-readout geometry)
    - merge_m4k_one_img






## kuiper_reduction package

In order to make the kuiper_reduce module more friendly with RTS2, we
created the kuiper_reduction package. This puts different parts 
of the analysis in different modules in case you don't want to import 
the entire analysis package. Also, its functions and methods operate
directly on astropy fits fds rather than filenames. 


```
from kuiper_reduction.merge import merge_fitsfd
merged_fd = merge_fitsfd( fitsfd )
```
