from rts2 import scriptcomm
from rts2 import sextractor
from scottSock import scottSock
sepPresent = False
try:
        import sep
        sepPresent = True
except Exception as ex:
        pass

import os
from pylab import *
from scipy import *
from scipy import optimize
import numpy
import pickle



def getTries(pathway="/home/rts2obs/rts2images"):
	focusfiles = [x for x in os.listdir(pathway) if "foc_20171204" in x and ".fits" in x]
	tries = {}
	for f in focusfiles:
		#print f
		num = f.split('_')[2].split(".")[0]
		tries[float(num)] = pathway+"/"+f

	return tries

def __sepFindFWHM(tries):
	from astropy.io import fits
	import math
	import traceback
	focpos=[]
	fwhm=[]
	fwhm_min=None
	fwhm_MinimumX=None
	keys = list(tries.keys())
	keys.sort()
	ln2=math.log(2)
	for k in keys:
		try:
			fwhms=[]
			ff=fits.open(tries[k])
			# loop on images..
			for i in range(1,len(ff)-1):
				data=ff[i].data
				bkg=sep.Background(numpy.array(data,numpy.float))
				sources=sep.extract(data-bkg, 5.0 * bkg.globalrms)
				#self.log
				print('I','bkg gobalrms {}'.format(bkg.globalrms))

				for s in sources:
					fwhms.append(2 * math.sqrt(ln2 * (s[15]**2 + s[16]**2)))


			im_fwhm=numpy.median(fwhms)
			# find median from fwhms measurements..
			
			#self.log
			print('I','median fwhm {}'.format(numpy.median(fwhms)))
			#self.log
			print('I','offset {0} fwhm {1} with {2} stars'.format(k,im_fwhm,len(fwhms)))
			focpos.append(k)
			fwhm.append(im_fwhm)
			if (fwhm_min is None or im_fwhm < fwhm_min):
				fwhm_MinimumX = k
				fwhm_min = im_fwhm
		except Exception as ex:
			#self.log
			print('W','offset {0}: {1} {2}'.format(k,ex,traceback.format_exc()))

	#self.log('I','pickling')
	#fd = open( "rts2.pkl", 'w' )
	#pickle.dump(sources, fd)
	#fd.close()
	return focpos,fwhm,fwhm_min,fwhm_MinimumX


def main():
	tries = getTries()
	focpos, fwhm, fwhm_min, fwhm_MinimumX = __sepFindFWHM(tries)
	print focpos, fwhm, fwhm_min, fwhm_MinimumX

main()
