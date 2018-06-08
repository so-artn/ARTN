#!/usr/bin/python
from astropy.io import fits
from rts2 import scriptcomm
from scottSock import scottSock

sepPresent = False
try:
	import sep
	sepPresent = True
except Exception as ex:
	pass


class ShiftFocus (scriptcomm.Rts2Comm):
	"""Take and process focussing data."""

	def __init__(self,exptime = 10,step=20,attempts=10,filterGalaxies=False):
		scriptcomm.Rts2Comm.__init__(self)
		self.log('I', 'This is a test')
		self.exptime = exptime
		self.step = step
		self.focuser = "F0"
		self.attempts = attempts

		# if |offset| is above this value, try linear fit
		self.linear_fit = self.step * self.attempts / 2.0
		# target FWHM for linear fit
		self.linear_fit_fwhm = 3.5
		self.filterGalaxies = filterGalaxies

	def shift_focus( self ):
		self.setValue("shiftfocus", 1, "C0")
		self.setValue( "exposure", 5 )
		img = self.exposure(self.before_readout, '%b/shiftfoc_%N_.fits' )
		#to_dataserver( img )
		

	def before_readout(self):
		
		self.log('I','Before readout')

	def run(self):
		self.shift_focus()

def to_dataserver( fname, outfile='test.fits', clobber=True ):

	fitsfd = fits.open( fname )
		

        width = 0
        height = 0
        for ext in fitsfd:
                if hasattr( ext, 'data' ):
                        if ext.data is not None:
                                width+=ext.data.shape[0]
                                height+=ext.data.shape[1]

        fitsfd.close()
        fsize = os.stat(fname).st_size

        fd = open(fname, 'rb')


        if clobber:
                clobber_char = '!'
        else:
                clobber_char = ''
        meta = "          {} {}{} 1 {} {} 0".format( fsize, clobber_char, '/home/rts2obs/shiftfocus.fits', width, height )
        meta = meta + (256-len(meta))*' '

        data = meta+fd.read()
        lendata = len(data)
        soc = scottSock( '10.30.1.1', 6543 )

        counter = 0
        socsize = 1024
        buffsize = 0
        while buffsize < len(data):
                sent = soc.send( data[buffsize:buffsize+1024] )
                buffsize+=sent


f=ShiftFocus()

f.run()

