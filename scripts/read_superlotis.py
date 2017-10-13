import os
import subprocess
import numpy

#uvot - 1a VBU 
#Stand - standard stars
#AzTec - VBR 60 60 60 x3

nami = 11
rai = 4
deci = 5
typi = 14
tid_index = 0 

filt_dict = {0:"120", 1:"60", 2:"60", 3:"60", 4:"60", 5:"60"}
type_dict = {"UVOT":0, "AzTEC":1}

class lotisimport:
	def __init__(self, sr):
		os = findoffset(sr)
		self.name = sr[nami-os]
		self.ra = formatcoord(sr[rai-os])
		self.dec = formatcoord(sr[deci-os])
		self.type = type_dict[sr[typi-os]]
		print "\nINIT: {}, {}, {}, {}".format(self.name, self.ra, self.dec, self.type)

def findoffset(sr):
	for key in type_dict.keys():
		if key in sr:
			typeindex = sr.index(key)
			return typi-typeindex
	return 0 

def formatcoord(coord):
	return "{}:{}:{}".format(coord[0:2],coord[2:4],coord[4:])

def findtargetid_targetlist(lotis_obj):
	cmd = "rts2-targetlist | grep -i {} > tmp.txt".format(lotis_obj.name)
	print cmd
	subprocess.call(cmd, shell=True)
	fi = open("tmp.txt")
	for row in fi:
		splitrow = row.split()
		#for ii,i in enumerate(splitrow):
		#	print ii, i
		if any(x for x in splitrow if str.lower(lotis_obj.name) in str.lower(x)):
			print "found index"
			targetid = splitrow[tid_index]
			return targetid
	fi.close()
	cmd = "rm tmp.txt"
	subprocess.call(cmd, shell=True)
	createobj(lotis_obj)

def setscript(obs_type, targetid):
	script = "BIG61.OFFS=(2m,0) "
	if obs_type == 0: #uvot
		for filt in [3,2,0]:
			tmp = "filter={} E {} ".format(str(filt), filt_dict[filt])
			script += (tmp * 3) 
	if obs_type == 1: #aztec
		for filt in [3,2,1]:
                        tmp = "filter={} E {} ".format(str(filt), filt_dict[filt])
                        script += (tmp * 3)
	cmd = "rts2-target -c C0 -s \"{}\" {}".format(script, targetid)
	print cmd
	subprocess.call(cmd, shell=True)

def createobj(lotis_obj):
	#THIS DOESN"T WORK NEED TO TALK TO SCOTT
	cmd = "python /home/scott/git-clones/rts2/scripts/newtarget.py --create {} {} {}".format(lotis_obj.ra, lotis_obj.dec, lotis_obj.name)
	print cmd
	subprocess.call(cmd, shell=True)
	return findtargetid_targetlist(lotis_obj)


def set_queue(targetids):
	if len(targetids) > 0:
		targetstring = ""
		for iid in targetids:
			targetstring += " {}".format(iid)
		cmd = "rts2-queue --queue plan --clear{}".format(targetstring)
		print cmd
		#subprocess.call(cmd, shell=True)
	
lotisdata = "lotis.txt"
fi = open(lotisdata, "r")
targetids = []
names = []
for row in fi:
	splitrow = row.split()
	if len(splitrow) > 14:
		#try:
		obj = lotisimport(splitrow)
		if obj.type in type_dict.values() and obj.name not in names:
			targetid = findtargetid_targetlist(obj)
			if targetid != 0: #remove once I figure out how to get newtarget working
				setscript(obj.type, targetid)
				targetids.append(targetid)
				names.append(obj.name)
		#except:
		#	print "bust"
		#	pass
setq_input = raw_input("Set Queue with LotisData? [y/n]: ")
if setq_input in ["Y","y","yes"]:
	set_queue(targetids)

print "fin"

	#call rts2-targetlist | grep -a name > file
	#read file 
	#test if an observation is in there
	#if yes: get targetid
	#	 set script automoatically to UBV:120, 60 60 x3

	#if not: call create script
	#then find targetid
	#then set script automatically to VBU:60, 60, 120 x3

	#live your life: love your family
	

