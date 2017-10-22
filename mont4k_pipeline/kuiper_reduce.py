#!python
#
# Test Kuiper 61" data reduction with python ccdproc
#
# Trying to follow along with the ccdproc example at
# http://nbviewer.jupyter.org/gist/mwcraig/06060d789cc298bbb08e
# with modifications to handle multi-extension FITS files
#
# I wrote code to try to group images into bias, flat,
# twiflat, object - it's in ~/python/ccdreduce/makelists.py
# and have inserted that into this file now.

# Benjamin Weiner, bjw@as.arizona.edu

# To run this, to process all the files in the current directory, try:
# python
# import sys 
# sys.path.append('/Users/bjw/stellarview/python')
# import kuiper_reduce as kp
# fstruct = kp.getheaders('')
# imagelist = []
# for image in fstruct:
#    imagelist.append(image['filename'])
# kp.process_and_merge_list(imagelist)


# Requires astropy, ccdproc, maybe msumastro
# maybe bottleneck

# needed for py notebook only?
# %matplotlib inline

import matplotlib.pyplot as plt
import numpy as np

from astropy.modeling import models
from astropy import units as u
from astropy import nddata
from astropy.io import fits

import ccdproc
import os

# For faster median of masked arrays until numpy 1.9
import bottleneck as bn

# May be no longer needed?

from msumastro import ImageFileCollection, TableTree

nddata.conf.warn_unsupported_correlated = False

# Adjust for instrument, or calculate values?
gain = 3.17 * u.electron / u.adu
readnoise = 6.4 * u.electron

#####
# convenience functions

# median of single array using bottleneck
def bn_median(masked_array, axis=None):
    data = masked_array.filled(fill_value=np.NaN)
    med = bn.nanmedian(data, axis=axis)
    # construct a masked array result, setting the mask from any NaN entries
    return np.ma.array(med, mask=np.isnan(med))

# median of a stack of masked arrays
def avg_over_images(masked_arr, axis=0):
    """
    Calculate average pixel value along specified axis
    """
    return ma.mean(masked_arr, axis=axis)

def med_over_images(masked_arr, axis=0):
    """
    Calculate median pixel value along specified axis    
    Uses bottleneck.nanmedian for speed
    """
    dat = masked_arr.data.copy()
    dat[masked_arr.mask] = np.NaN
    return bn.nanmedian(dat, axis=axis)

# convenient plotting.  Set min max stretch to scale*stddev
def plot_image(image, scale):
    min, max, mean, stddev = imstats(np.asarray(image))
    plt.figure(figsize=(15, 15))
    plotmax = mean + scale*stddev
    plotmin = mean - scale*stddev
    plt.imshow(image, vmax=plotmax, vmin=plotmin)
    # need to write the plot out somehow
    return
    
# image statistics
imstats = lambda dat: (dat.min(), dat.max(), dat.mean(), dat.std())

# read all data in a directory
def read_data(direc):
    data_dir = direc
    images = ImageFileCollection(data_dir, keywords='*')
    return images

# write an image or a list of images. prefix can be eg 'r_' or 'reduced_' or nothing
# Argument is a list of image objects not filenames?
def save_images(images,prefix):
    for im in images:
        im_hdu_list = im.to_hdu()
        fname_base = os.path.basename(im.header('filename'))
        newname = prefix + fname_base
        try:
            im_hdu_list.writeto(newname)
        except: 
            mylog('Failed write, trying to overwrite: '+newname)
            im_hdu_list.writeto(newname,overwrite=True)
    return

# Convenience function for logging.  Could be print or output to file
def mylog(mystring):
    print mystring
    return

#####
# Get image lists / info


# get information from image headers - extract common header fields
# and return a table or similar

# fields to get:
# filename, title, nx, ny, nextend,
# ut, ra, dec, epoch, alt, az, airmass, filter, exptime,
# x bin, y bin, imtype,
# mean, rms

# import fits

#from astropy.table import Table
#import astropy.io.fits as fits

# if filelist='' or allfitsfiles is True then find all *.fits files in
# directory
def getheaders(filelist,dir='',allfitsfiles=False):

    if dir=='':
        dir = os.getcwd()
    if allfitsfiles==True or allfitsfiles>0 or filelist=='':
        allfilelist = os.listdir(dir)
        filelist=[]
        for name in allfilelist:
            if os.path.isfile(name) and name.endswith('.fits'):
                filelist.append(name)
    
    # read in header translation fields as dictionary
    # initialize dictionary structure
    outstruct = []
    nfiles = len(filelist)
    
    for file in filelist:
        fullname = os.path.join(dir,file)
        # use extension 0 or 1? 0 for header, 1 for data?
        try:
            hdus = fits.open(fullname)
            # get file header
            hdr = hdus[0].header
            # if there is more than 1 hdu, use hdu=1 for the data
            ndhus = len(hdus)
            # ihdudata = 0
            #if nhdus > 1 :
            ihdudata = 1
                 # hdr1 = hdus[ihdudata].header
            # do simple image statistics to get mean and rms of image
            # No data may yield a nan
            mean1 = np.mean(hdus[ihdudata].data)
            # median1 = np.median(hdus[1].data)
            rms1 = np.std(hdus[ihdudata].data)        
            # get each keyword and put into structure
            struct1 = {'filename':file, 'title':hdr['object'], 'nextend':hdr['nextend'],
                    'ut':hdr['ut'], 'ra':hdr['ra'], 'dec':hdr['dec'],
                   'equinox':hdr['equinox'], 'azimuth':hdr['azimuth'], 'elev':hdr['elevat'],
                   'airmass':hdr['airmass'], 'filter':hdr['filter'], 'exptime':hdr['exptime'],
                   'imtype':hdr['imagetyp'], 'mean':mean1, 'rms':rms1 }

                    #'nx':hdr['naxis1'],'ny':hdr['naxis2'], 'bin':hdr['ccdsum'],
            outstruct.append(struct1)
            hdus.close()
        except:
            print 'Error opening file: ',fullname
        
    # more stuff?
    mylog("getheaders: Got headers for {0} files".format(nfiles))

    return outstruct


# Given a structure of common header keywords and image mean and
# rms, try to figure out which images are which and make lists
# of their file names.

def makelists(imstruct):

    zerolist = []
    darklist = []
    domelist = []
    skylist = []
    objectlist = []
    badlist = []
    stdlist = []

    zeromax = 3000.0
    darkmax = 5000.0
    flatmin = 5000.0
    objectmax = 5000.0
    skyflatmin = 8000.0

    # read in optional parameters for instrument to override defaults
    # for flatmin, etc?
    
    for image in imstruct:
        imtype = image['imtype']
        fname = image['filename']
        title = image['title']
        exptime = image['exptime']
        immean = image['mean']
        if (imtype == 'zero' or imtype == 'bias') and exptime<0.001:
            if immean < zeromax:
                zerolist.append(fname)
            else:
                print 'makelists Warning: zero with counts>',zeromax,fname
                badlist.append(fname)
        elif imtype == 'dark' and exptime>0.1:
            if immean < darkmax:
                darklist.append(fname)
            else:
                print 'makelists Warning: dark with counts>',darkmax,fname
                badlist.append(fname)
        elif imtype == 'flat' and exptime>0.1:
            if immean > flatmin:
                domelist.append(fname)
            else:
                print 'makelists Warning: flat with counts<',flatmin,fname
                badlist.append(fname)
        elif imtype == 'object':
            if immean < objectmax:
                objectlist.append(fname)
                # check for Landolt or standard in image title
                ifstd = ('landolt' in title.lower() or 'standard' in title.lower())
                if ifstd:
                    print 'makelists: standard field ',fname
                    stdlist.append(fname)
            else:
                print 'makelists Warning: object with counts>',objectmax,fname
                if immean > skyflatmin:
                    skylist.append(fname)
                    print 'makelists: using as skyflat object with counts>',skyflatmin,fname
                else:
                    badlist.append(fname)
        else:
            print 'makelists: Warning: didnt understand image type for ',fname
            badlist.append(fname)

    # print lists to files?
    write_list_to_file(zerolist,'zero.list')
    write_list_to_file(domelist,'dflat.list')
    write_list_to_file(skylist,'sflat.list')
    write_list_to_file(objectlist,'object.list')
    write_list_to_file(stdlist,'standard.list')
    write_list_to_file(badlist,'badimages.list')

    mylog("makelists: Made the *.list files")

    # subset lists by binning and filters?
    
    return zerolist,domelist,skylist,objectlist,stdlist,badlist

# This will create/overwrite files
def write_list_to_file(filelist,fname):
    f = open(fname,'w')
    for name1 in filelist:
        line = '%s\n' % (name1)
        f.write(line)
    f.close()
    return

#####
# Processing functions

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
     oscan1 = ccdproc.subtract_overscan(data1, fits_section=data1.header['BIASSEC'], add_keyword={'overscan': True, 'calstat': 'O'}, model=models.Polynomial1D(1))
     trim1 = ccdproc.trim_image(oscan1, fits_section=oscan1.header['TRIMSEC'], add_keyword={'trimmed': True, 'calstat': 'OT'})
     fits.update(fname, trim1.data, header=trim1.header, ext=i)
  hdulist.close()
  mylog("Overscan and trim {0}".format(fname))
  return

# overscan and trim a list
def oscan_trim_list(filelist):
    for fname in filelist:
        oscan_trim_file(fname)
    return
        
# overscan, trim, combine with sigma clipping

# combine a multi ext fits file. First we get all the data into combine_list,
# then rearrange it so that all the 1st extensions are combined, all the 2nd,
# and so on.
# can we pass in an option to scale during combine?
def combine_list_to_file(listname,outname,read_from_file=False,combine='median'):
    if read_from_file == True:
        listfile = open(listname,'r')
        flist = listfile.read().splitlines()
        listfile.close()
    else:
        flist = listname
    if len(flist) == 0:
        mylog('combine_list got a list of zero files '+listname)
        return
    # combine_list will be a list of lists of data objects
    combine_list = []
    ifirst = True
    for line in flist:
    	fname = line.strip()
	hdulist = fits.open(fname)
	if ifirst == True:
            try:
                hdulist.writeto(outname)
            except:
                mylog('Failed write, trying to overwrite: '+outname)
                hdulist.writeto(outname,overwrite=True)
            ifirst = False
        nhdus = len(hdulist)
        # Data of this image will be listed in a single entry in combine_list
	tmplist = []
        if nhdus > 1:
            istart = 1
        else:
            istart = 0
	for i in range(nhdus)[istart:] :
            data1 = ccdproc.CCDData(hdulist[i].data, unit="adu")
            head1 = hdulist[i].header
            head1['filename'] = fname
            data1.header = head1
            tmplist.append(data1)
            # combine_list.append(ccdproc.CCDData(data=data1, meta=head1, unit="adu"))
        hdulist.close()
      	combine_list.append(tmplist)
    for	i in range(nhdus)[istart:] :
        # Take the i-1'th data extension (0-based) from each image and
        # put these into a list to combine
        tmplist = [ elem[i-1] for elem in combine_list ]
	combo = ccdproc.Combiner(tmplist)
        if combine == 'average':
            output1 = combo.average_combine()
        else:
            output1 = combo.median_combine()
        # Header output by combine is minimal and causes failure
        # fits.update(outname, output1.data, header=output1.header, ext=i)
        # but outputting without header field leaves only 6 line header
        # fits.update(outname, output1.data, ext=i)
        # copying header from the first input image works.
        fits.update(outname, output1.data, header=tmplist[0].header, ext=i)
    mylog("Combined to output file {0}".format(outname))
    return

# bias correct a file with a master bias file
def bias_file_by_file(fname,biasname):
    hdulist = fits.open(fname)
    hdubias = fits.open(biasname)
    nhdus = len(hdulist)
    if nhdus > 1:
        istart = 1
    else:
        istart = 0
    for i in range(nhdus)[istart:] :
        data1 = ccdproc.CCDData(hdulist[i].data, unit="adu")
        data1.header = hdulist[i].header
        bias1 = ccdproc.CCDData(hdubias[i].data, unit="adu")
        bias1.header = hdubias[i].header
        commentstr = "Bias image is "+biasname
        # proc1 = ccdproc.subtract_bias(data1,bias1,add_keyword={'bias': True, 'calstat': 'OTZ', 'history':commentstr} )
        proc1 = ccdproc.subtract_bias(data1,bias1,add_keyword={'bias': True, 'calstat': 'OTZ'} )
        fits.update(fname, proc1.data, header=proc1.header, ext=i)
        # fits.update(fname, proc1.data, ext=i)
    hdulist.close()
    hdubias.close()
    mylog("Bias corrected {0} with {1}".format(fname,biasname))
    return

# apply bias subtraction to a list
def bias_list_by_file(filelist, biasname):
    for fname in filelist:
        bias_file_by_file(fname, biasname)
    return

# flat correct a file with a master flat file
# This does not try to fix the flat normalization to a common value
# across the extensions, because that isn't implemented in flat_correct yet
def flat_file_by_file(fname,flatname):
    hdulist = fits.open(fname)
    hduflat = fits.open(flatname)
    nhdus = len(hdulist)
    if nhdus > 1:
        istart = 1
    else:
        istart = 0
    for i in range(nhdus)[istart:] :
        data1 = ccdproc.CCDData(hdulist[i].data, unit="adu")
        data1.header = hdulist[i].header
        flat1 = ccdproc.CCDData(hduflat[i].data, unit="adu")
        flat1.header = hduflat[i].header
        if i==1:
            flatscale = np.mean(flat1)
        # flat1 = flat1/flatscale
        commentstr = "Flat image is "+flatname+" with scale"+str(flatscale)
        proc1 = ccdproc.flat_correct(data1,flat1,add_keyword={'flat': True, 'calstat': 'OTZF', 'history':commentstr} )
        fits.update(fname, proc1.data, header=proc1.header, ext=i)
        # fits.update(fname, proc1.data, ext=i)
    hdulist.close()
    hduflat.close()
    mylog("Flat corrected {0} with {1}".format(fname,flatname))
    return

# combine biases into a master bias
# using combine_list_to_file allows a multiple fits extension image
def make_master_bias(imagelist,bias_name='Bias_master.fits'):
    combine_list_to_file(imagelist,bias_name,combine='median')
    return

# combine one list of flat names into a master flat
# It would be useful to have an option to scale the flats during combine
# Combining flats with average lessens the need to scale, but
# scaling and combining skyflats with median is desirable to reject stars.
def make_master_oneflat(imagelist,flat_name='Flat_master.fits',combine='average'):
    if combine == 'median':
        combine_list_to_file(imagelist,flat_name,combine='median')
    else:
        combine_list_to_file(imagelist,flat_name,combine='average')
    return

# take a list of flats, sort it by filter, and combine each set into
# a master flat for that filter.
# This assumes that all HDUs of an image are in the same filter. If you
# have a multicolor camera where each HDU is a diff filter but they
# are always fixed, can use make_master_oneflat directly.
# This could have used the structure made by getheaders, but it doesn't.
# Return a list of the master flat names
def make_master_flats(listname,flat_prefix='Flat_',filter_key='FILTER',read_from_file=False):
    if read_from_file == True:
        listfile = open(listname,'r')
        flist = listfile.read().splitlines()
        listfile.close()
    else:
        flist = listname
    filter_name_dict, filters_uniq = find_filters_list(flist,filter_key)
    mylog("Found flats in filters {0}".format(filters_uniq))
    list_outnames = []
    
    for name in filters_uniq:
        # find matching filenames and make list
        # flatlist1 = filter_name_dict[name]
        flatlist1 = [fname for fname, filt in filter_name_dict.items() if filt == name]
        # remove any spaces from filter name
        filtnamesquash = name.replace(" ","")
        outname = flat_prefix + filtnamesquash + '.fits'
        make_master_oneflat(flatlist1, flat_name=outname)
        list_outnames.append(outname)
        mylog("Combined master flat {0}".format(outname))
    # all done
    return list_outnames

# Find the filter names of each image and return a dictionary
# of image: filter_name, and a list of unique filter names
def find_filters_list(flist,filter_key='FILTER'):
    filter_name_dict = {}
    filter_names = []
    
    for line in flist:
    	fname = line.strip()
	hdulist = fits.open(fname)
        filtname = hdulist[0].header[filter_key]
        if filtname == '' and len(hdulist)>1 :
            filtname = hdulist[1].header[filter_key]
        filter_names.append(filtname)
        #struct1 = {'filename': fname, 'filter': filtname}
        # struct1 = {filtname: fname}
        # struct1 = {fname: filtname}
        filter_name_dict[fname] = filtname
        hdulist.close()
    # uniq-ify the filter name list
    filters_uniq = list( set(filter_names) )
    # this might also give a unique set 
    filters_uniq = list( set( filter_name_dict.values() ) )
    mylog("find_filters_list made a dictionary and found filters: ")
    # print filter_name_dict
    print filters_uniq
    return filter_name_dict, filters_uniq

# ccdproc.flat_correct automatically normalizes by the mean of the flat,
# ie it divides by the flat and multiplies back by the mean.
# When proc'ing multiple extensions separately, we really want them all
# to be normalized by the same factor. So we need to calculate these
# flat means by extension and divide each image ext by its mean, then
# multiply up by the mean of means.
def renormalize_by_flat(image,flat,read_from_file=False):
    if read_from_file == True:
        hduimage = fits.open(image)
        hduflat = fits.open(flat)
    else:
        hduimage = image
        hduflat = flat
    nhdus = len(hdulist)
    # Do nothing if image only has 0-1 data extension
    if nhdus <=1:
        mylog("Don't need to renormalize a 1-extension flat, returning")
        return
    secmeans = np.zeros(nhdus-1)
    for i in range(nhdus)[1:] :
        secmeans[i-1] = hduflat[i].data.mean()
    totmean = secmeans.mean()
    secmeans = secmeans / totmean
    mylog("Renormalized flat means by extension: {0} ".format(secmeans))
    for i in range(nhdus)[1:] :
        #hduimage[i].data = hduimage[i].data / secmeans[i-1]
        hduimage[i].divide(secmeans[i-1] * hduimage[i].unit)
        if read_from_file == True:
            fits.update(image, hduimage[i].data, header=hduimage[i].header, ext=i)
    if read_from_file == True:
        hduimage.close()
        hduflat.close()
    else:
        # This may be unnecessary if these are already pointing at identical structure
        image = hduimage
    # done?
    return
    
# make flats w/o dark subtraction
# these are combining dome flats most likely
# arguments are the collection of images, and a filter string like 'R'
# def make_flat_filter(images, filtername):

# Do something with twilight flats?

# Construct a supersky from night images, or a fringe frame?

# Process science images
# master_bias and master_flat are input params so that you can grab
# them from a different directory, etc
# assumes your objects are all in same filter as the flat
#
def proc_objects_filter(imagelist, master_bias, master_flat, filtername=''):
    # list of objects?
    # pick only those matching filtername?
    for fname in imagelist:
        bias_file_by_file(fname, master_bias)
        flat_file_by_file(fname, master_flat)
    mylog("Processed image {0}".format(fname))
    return

# Take a list of objects and a list of flats and sort so that
# objects get the correct flat
def proc_objects_all(imagelist, master_bias, flatlist):
    # sort by matching objects and flats?
    filter_name_dict, filters_uniq = find_filters_list(imagelist)
    flat_name_dict, flat_filters_uniq = find_filters_list(flatlist)
    for filter1 in filters_uniq:
        objlist1 = [fname for fname, filt in filter_name_dict.items() if filt == filter1]
        if filter1 in flat_filters_uniq:
            flats_match = [fname for fname, filt in flat_name_dict.items() if filt == filter1]
            # if more than one flat matches, which is unlikely, take the first
            flat_match = flats_match[0]
            proc_objects_filter(objlist1, master_bias, flat_match)
            mylog("processed objects for {0}".format(filter1))
        else:
            mylog("proc_objects_all: couldnt find flat for {0}".format(filter1))
    # end
    return

# given a list of images, divide them up into biases, flats, objects,
# make master biases and flats, then apply them
def process_list(imagelist):

    imstruct = getheaders(imagelist)
    zerolist,domelist,skylist,objectlist,stdlist,badlist = makelists(imstruct)

    masterbiasname = 'Bias_master.fits'
    oscan_trim_list(zerolist)
    make_master_bias(zerolist,bias_name=masterbiasname)

    oscan_trim_list(domelist)
    bias_list_by_file(domelist,masterbiasname)
    flatnames = make_master_flats(domelist,flat_prefix='Flat_',filter_key='FILTER')
    # Are we going to use domeflats or skyflats?
    oscan_trim_list(skylist)
    bias_list_by_file(skylist,masterbiasname)
    skyflatnames = make_master_flats(skylist,flat_prefix='Skyflat_',filter_key='FILTER')

    oscan_trim_list(objectlist)
    # flatfield with domes
    proc_objects_all(objectlist, masterbiasname, flatnames )
    # flatfield with skies only
    # proc_objects_all(objectlist, masterbiasname, skyflatnames )
    mylog("Processed images from list {0}".format(imagelist))
    return objectlist

# merge the two extensions of a list of mont4k image files
# save the raw images?
# Some of this is taken from Scott Swindell's m4kproc.py
def merge_m4k_list(imagelist,raw_dir="raw",merged_dir="merged"):
    # make the raw or the merged directory?
    if merged_dir.endswith("/"):
        merged_dir = merged_dir[:-1]
    if not os.path.exists(merged_dir):
        os.mkdir(merged_dir)
    for fname in imagelist:
        hdulist = fits.open(fname)
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
        fname_base = os.path.basename(fname)
        newname = merged_dir + "/" + fname_base
        # newname = "{0}/{1}".format(merged_dir, fname_base)
        try:
            hdulist_new.writeto( newname )
        except:
            mylog('Failed write, trying to overwrite: '+newname)
            hdulist_new.writeto( newname, overwrite=True )
        hdulist_new.close()
        mylog("merge_m4k_list: wrote merged image {0}".format(newname))
    #
    return

# Process all the images and merge only the objects
def process_and_merge_list(imagelist):
    objectlist = process_list(imagelist)
    merge_m4k_list(objectlist)
    return objectlist    

# bad pixel mask

# cosmic ray cleaning

# make plots
