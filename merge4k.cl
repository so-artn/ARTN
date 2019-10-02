# This is a Sep 2015 version on gerard to convert Mont4k AZCam extended fits
# images to normal fits images that all IRAF tasks can deal with
#
# It should work for images taken prior to Sept 2015 (?)
#
# Inputs an individual raw mont4k extended fits (file binned 1x1, 
# 2x2, 3x3, or 4x4) # and outputs a merged normal fits image
#
# A previously merged image is left unchanged by this routine
# 
# The input images are saved unchanged in a subdirectory named raw
# (The exact image size is read from the headers (so this routine
# works for different numbers of rows and columns -- as long as the
# two sides of the ccd have the same number of columns!)
# !!  Does not yet work on any images with detector rows not from 1 to 4096
# development!)
#
# The only processing that has been done on the output image is the
# correction for electronic cross-talk, the overscan subtraction and 
# trimming of each side required to merge the two halves into one image  

procedure merge4k (rawimg)

string rawimg      {prompt="raw multi-extension mont4k image(s)"}
struct  *imagelist

begin
     string rimgstr,rimg,pkg,rimg0,rimg1,rimg2,img1,img2,imagel
     string hdr,typ,bin,ysec0,ysec,biassec,datasec,datasec2,datasec3
     string xt1,xt2,ccdsec,ccdsec1,strval
     real xtcoef,ovsc1,ovsc2
     int i,j,n1,n2,nlast,x1,x2,ibin,intval

     cache imgets

# 2005 original cross-talk coefficient
#    xtcoef = 0.00264
# Aug 2008 cross-talk coefficient, after pre-amps installed
#     xtcoef = 0.00026
# Sep 2014 cross-talk coefficient
     xtcoef = 0.00080


# check to make sure fitsutil package is loaded; needed to 
#  split extended fits file with fxsplit
     if (access("pkg")) del ("pkg")
     lpar ("fxextract", >& "pkg")
     head ("pkg", nlines=1) | scan (pkg)
     if (substr(pkg,1,5) == "ERROR")  \
       error(1,"Load fitsutil package in order to run merge4k")
     del ("pkg")

# check to make sure ctio package is loaded; needed to create 
#  new normal fits file with original headers using imcreate
     lpar ("imcreate", >& "pkg")
     head ("pkg", nlines=1) | scan (pkg)
     if (substr(pkg,1,5) == "ERROR")  \
       error(1,"Load ctio package in order to run merge4k")
     del ("pkg")

# create the raw directory if it doesn't exist
     if (!access("raw")) mkdir ("raw")

     nhedit.update = yes
     nhedit.add = yes
     nhedit.addonly = no
     nhedit.rename = no
     nhedit.verify = no
     nhedit.show = no

     colbias.median = yes
     colbias.interactive = no
     
     rimgstr = rawimg
     if (substr(rimgstr,1,1) != "@") {
       i = strlen(rimgstr)
       if (substr(rimgstr,i-3,i) != "fits") rimgstr = rimgstr//".fits"
     }
     imagel = mktemp ("img")
     files (rimgstr, > imagel)

# read and process each image in the list
     imagelist = imagel
     while (fscan (imagelist,rimg) != EOF) {

# only try to convert an image if it exists, hasn't already been converted
# (ie image type is not real), and the data type is ushort (ie it probably
# came direct from AZCam)
     if (access (rimg)) {
       hdr=""
       imhead ((rimg//"[0]")) | scan (hdr)
       i = strlen(hdr)
       typ = substr(hdr,i-5,i-2)
       if(typ != "real") {
         imhead ((rimg//"[1]")) | scan (hdr)
         i = strlen(hdr)
         typ = substr(hdr,i-7,i-2)
       }
       if(typ == "ushort") {
         i = strlen(rimg)
         rimg0 = substr(rimg,1,i-5)//"0"//".fits"
         rimg1 = substr(rimg,1,i-5)//"1"//".fits"
         rimg2 = substr(rimg,1,i-5)//"2"//".fits"

# delete any images with this same name in the raw directory
         if (access ("raw/"//rimg)) del (("raw/"//rimg), verify=no)
# delete any previously existing image parts & temp files in the current directory
         if (access (rimg0)) del (rimg0,verify=no)
         if (access (rimg1)) del (rimg1,verify=no)
         if (access (rimg2)) del (rimg2,verify=no)
         if (access ("stupid")) del ("stupid", verify=no)

# split the extended fits input image into its parts
         fxsplit (rimg, verb=no)

# do the crosstalk correction for each half
         xt1 = mktemp ("xt1")
         xt2 = mktemp ("xt2")
         imar ((rimg//"[1]"), "*", xtcoef, xt2)
         imar ((rimg//"[2]"), "*", xtcoef, xt1)
         imar (rimg1, "-", xt1, rimg1)
         imar (rimg2, "-", xt2, rimg2)
         imdel (xt1,verify=no)
         imdel (xt2,verify=no)

# get binning and image dimension info
         imgets ((rimg//"[0]"), "ccdsum")
         bin = imgets.value
         ibin = int(substr(bin,1,1))

         imgets ((rimg//"[1]"), "ccdsec")
         ccdsec = imgets.value
         j = stridx(",",ccdsec)
         ysec0 = substr(ccdsec,j,strlen(ccdsec))

# construct various image section sizes
         imgets ((rimg//"[1]"), "datasec")
         datasec = imgets.value
         i = stridx(":",datasec)
         j = stridx(",",datasec)
         n1 = int(substr(datasec,i+1,j-1))
         ysec = substr(datasec,j,strlen(datasec))

         hselect ((rimg//"[1]"), "naxis[1]", yes) | scan (nlast)
         x1 = n1+6; x2 = nlast
         biassec = "["//str(x1)//":"//str(x2)//ysec
         x1 = 1; x2 = n1
         datasec2 = "["//str(x2)//":"//str(x1)//ysec
         x1 = n1+1; x2 = 2*n1
         if (ibin == 3) {
           x1 = x1+1
           x2 = x2+1
         }
         datasec3 = "["//str(x1)//":"//str(x2)//ysec

# ccdsec no longer useful for working with binned images after 8/28/10 because it 
#   doesn't take binning into account; write ccdsec1 keyword in header also
         x1 = 1 ; x2 = 4096
         ccdsec = "["//str(x1)//":"//str(x2)//ysec0
         i = 4096/ibin
         hselect (rimg1, "naxis[2]", yes) | scan (j)
         x1 = 1 ; x2 = i ; if (i == 3) x2 = i+1
         ccdsec1 = "["//str(x1)//":"//str(x2)//ysec

# overscan correct and trim each half
         colbias (rimg1,rimg1, bias=biassec, trim=datasec, func="spline3", order=2)
         colbias (rimg2,rimg2, bias=biassec, trim=datasec, func="spline3", order=2)

         imstat ((rimg//"[1]"//biassec),field="midpt",format=no) | scan (ovsc1)
         imstat ((rimg//"[2]"//biassec),field="midpt",format=no) | scan (ovsc2)

# save the original extended fits image in the raw subdirectory
         rename (rimg, ("raw/"//rimg))

# create new image with header of original [0]
         imcreate (rimg, naxis=2, naxis1=i, naxis2=j, header="copy", pixtype="real", \
            reference=rimg0)

# copy the crosstalk-corrected, overscan-corrected, trimmed images into a combined image
# the following lines are required to make iraf realize what size rimg1 is
         print (("imstat "//rimg1), > "stupid")
         cl < "stupid", >& "dev$null"
         del ("stupid")

         imcopy ((rimg1//datasec), (rimg//datasec), verbose=no)
         img1 = rimg2//datasec2 ; img2 = rimg//datasec3
         imcopy ((rimg2//datasec2), (rimg//datasec3), verbose=no)

# delete header keywords that are duplicates or no longer needed
         hedit (rimg, "detsize", delete=yes, verify=no, show=no, update=yes)
         hedit (rimg, "nextend", delete=yes, verify=no, show=no, update=yes)

# add useful header keywords from [1]/[2] that are not in [0] and/or 
#   fix previously existing header keywords in combined image
         nhedit (rimg, "extend", "F", "")
         imgets (rimg1, "bunit")
         strval = imgets.value
         nhedit (rimg, "bunit", strval, "Physical unit of array values")
         imgets (rimg1, "ccdname")
         strval = imgets.value
         nhedit (rimg, "ccdname", strval, "CCD name")
         imgets (rimg1, "detsize")
         strval = imgets.value
         nhedit (rimg, "detsize", strval, "Detector size")
         imgets (rimg1, "ccdsize")
         strval = imgets.value
         nhedit (rimg, "ccdsize", strval, "CCD size")
         imgets (rimg1, "ampsec")
         strval = imgets.value
         nhedit (rimg, "ampsec", strval, "Amplifier section")
         imgets (rimg1, "detsec")
         strval = imgets.value
         nhedit (rimg, "detsec", strval, "Detector section")
         imgets (rimg1, "biassec")
         strval = imgets.value
         nhedit (rimg, "biassec", strval, "Bias section")
         imgets (rimg1, "trimsec")
         strval = imgets.value
         nhedit (rimg, "trimsec", strval, "Trim section")
# datasec keyword from rimg1 or rimg2 is only half of the image; use ccdsec1
#   value in datasec so that ds9 command will display the whole binned image
#         imgets (rimg1, "datasec")
#         strval = imgets.value
         nhedit (rimg, "datasec", ccdsec1, "Data section")
         nhedit (rimg, "ccdsec", ccdsec, "CCD section")
         nhedit (rimg, "ccdsec1", ccdsec1, "CCD section with binning")
         imgets (rimg1, "ovrscan1")
         intval = int(imgets.value)
         nhedit (rimg, "ovrscan1", intval, "Overscan on axis 1")
         imgets (rimg1, "ovrscan2")
         intval = int(imgets.value)
         nhedit (rimg, "ovrscan2", intval, "Overscan on axis 2")
         imgets (rimg1, "prescan1")
         intval = int(imgets.value)
         nhedit (rimg, "prescan1", intval, "Underscan on axis 1")
         imgets (rimg1, "prescan2")
         intval = int(imgets.value)
         nhedit (rimg, "prescan2", intval, "Underscan on axis 2")
         imgets (rimg1, "equinox")
         intval = int(imgets.value)
         nhedit (rimg, "equinox", intval, "Equinox of WCS")
         imgets (rimg1, "wcsdim")
         intval = int(imgets.value)
         nhedit (rimg, "wcsdim", intval, "WCS Dimensionality")
         imgets (rimg1, "ctype1")
         strval = imgets.value
         nhedit (rimg, "ctype1", strval, "Coordinate type")
         imgets (rimg1, "ctype2")
         strval = imgets.value
         nhedit (rimg, "ctype2", strval, "Coordinate type")

         nhedit (rimg, "xtkcoef", xtcoef, "Crosstalk coefficient")
         nhedit (rimg, "overscan", biassec, "Overscan region used")
         nhedit (rimg, "ovscmed1", ovsc1, "Overscan median counts, side 1")
         nhedit (rimg, "ovscmed2", ovsc2, "Overscan median counts, side 2")

# delete the intermediate split images
         imdel ((rimg0),verify=no)
         imdel ((rimg1),verify=no)
         imdel ((rimg2),verify=no)
       }
      } else {
       print ((rimg//" doesn't exist in this directory"))
     }
     }

     delete (imagel, verify=no)
     flpr
     flpr
     print (".")

end
# really the end

