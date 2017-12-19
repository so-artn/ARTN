#!python
#
#
#  Prompt user or read a text file for info to make a new target
#  and write a shell script that executes RTS2 commands to configure
#  an RTS2 observation of the target.
#
#  No sanity checking on values is performed, so if you enter an invalid
#  RA, Dec, filter number, etc, it won't tell you.
#
#  The format of the text file is:
#  name ra dec queue_id_num filtnum1 exptime1 n_exp1 filtnum2 exptime2 n_exp2 ...
#  for example:
#  sn2017glx 19:43:40.3 56:06:36.3 1101  1 180 3  3 120 3  2 240 3
#
#  BJW, Nov 29-30 2017
#
#  DONE: read input from file rather than interactive
#  TODO: have it execute the script rather than asking the user to run w/sh
#  TODO: have rts2-newtarget use command line args
#  TODO: read the filter names from a file and translate name-> number
#  TODO: dithering?

import os
import sys
import rts2
# from astropy.io import ascii

class filter_set:
    """Class to simplify look up of filter by name and number
    

    uses the python [] operator to lookup filter number
    or name. If you give it the name it will return the 
    number and vice versa. it also uses aliases for the
    lookup. RTS2 and the Galil like to use long names
    like "Harris-U" observers like short names like "U"
    either format is acceptable as long as the alias 
    is in the dictionary below. 
    """


    # Filter Name Aliases. 
    # The keyword is the name of the filter as told by 
    # the galil and RTS2, the value in the dict is a tupple
    # of possible aliases for each filter
    alias = {
            "Harris-U": ("U"), 
            "Harris-R": ("R"), 
            "Harris-V": ("V"), 
            "Arizona-I": ("I"), 
            "Harris-"U"": ("B")
            "Schott-8612": ("Schott")  }


    def __init__(self, filters = None):
        
        if filters is None:
            self._filter_list = []

        elif type(filters) == list:
            self._filter_list = filters

        elif type(filters) == dict:
            # this assumes that the keywords of the dictionary are 
            # the fitler names and the value is the filter number. 

            
            #sort by filter number and reverse look up. 
            for key, value in sorted(filters.iteritems(), key=lambda (k,v): (v,k)):
                self._filter_list.append( key )

        elif type(filters) == str or type(filters) == unicode:
            self._filter_list = str(filters).split()

        else:
            raise TypeError("Unexpected filter type {}, type must be string, unicode, list or dict".format(type(filters)))


    def check_alias( self, alias ):
        
        for name, aliases in self.alias.iteritems():
            if alias.lower() == name.lower():
                return alias

            else:
                for al in aliases:
                    if al.lower() == alias.lower():
                        return name

        # we didn't find the alias
        return None
            
    def __str__(self):
        return "<filter_set: "+str(self._filter_list)+">"

    def __repr__(self):
        return self.__str__()

    def __getitem__(self, key):
        if type(key) == int:
            return self._filter_list[key]

        elif type(key) == str or type(key) == unicode:
            realname = self.check_alias(key)
            if realname is not None:
                return self._filter_list.index(realname)
        raise ValueError( "cannot find filter {}".format(key) )

        



                


            


# This sets a default filter order, but it may be different between runs,
# so be careful.  Current default order is U,R,B,V,I,Schott.
# It returns a dictionary with key=filter name string, value = filter number
# In the future we should read the order from a file.
def set_filter_dict():
    filter_dict = {}
    filter_dict['U'] = 0
    filter_dict['B'] = 2
    filter_dict['V'] = 3
    filter_dict['R'] = 1
    filter_dict['I'] = 4
    filter_dict['Schott'] = 5
    # Also put the numbers into the dictionary so specifying a filter by
    # number-string is legal, ie you can use either 'R' or '1' to mean filter 1.
    for i in range(6):
        filter_str = '%1i' % (i)
        filter_dict[filter_str] = i
    return filter_dict

# This reads the filter order from the RTS2 proxy
# and returns a dictionary similar to what 
# set_filter_dict does:

def get_filter_dict(prx_name, prx_passwd):
    proxy = rts2.rtsapi.createProxy( "http://bigpop:8889", prx_name, prx_passwd  )
    filters_str = proxy.getValue("W0", "filter_names", True)
    return filter_set( filters_str )

# Write the rts2-newtarget script. This is run taking input rather than
# command line arguments, so I use the shell "<< EOF" syntax.
def make_newtarget_prompt(fname):
    radec_str = raw_input('Enter RA and Dec as hh:mm:ss dd:mm:ss, separated by a space: ')
    name_str = raw_input('Enter object name, eg Landolt_92_249: ')
    queue_id_str = raw_input('Enter numeric queue id, eg 1200: ')
    f = open(fname,'a')
    f.write('rts2-newtarget << EOF\n')
    f.write(radec_str+'\n')
    f.write('s\n')
    f.write(queue_id_str+'\n')
    f.write(name_str+'\n')
    f.write('EOF\n')
    f.close()
    print "New target script written to ",fname
    queue_id_num = int(queue_id_str)
    return queue_id_num

# Write the rts2-newtarget script using args passed in rather than
# interactive. Rather than writing to a file, pass back a string
# that has newlines in it.
def make_newtarget_args(name,ra,dec,queue_id_num):
    cmd = 'rts2-newtarget << EOF\n'
    radec_line = '%s %s\n' % (ra,dec)
    id_line = '%s\n' % (queue_id_num)
    name_line = '%s\n' % (name)
    cmd = cmd + radec_line
    cmd = cmd + 's\n'
    cmd = cmd + id_line
    cmd = cmd + name_line
    cmd = cmd + 'EOF\n'
    return cmd  

def set_target_defaults():
    lunardist = '20:'
    airmass = ':2.2'
    camera = 'C0'
    offset = 'BIG61.OFFS=(1m,0)'
    cmd = 'rts2-target'
    cmd = cmd + ' -c ' + camera
    cmd = cmd + ' --lunarDistance ' + lunardist
    cmd = cmd + ' --airmass ' + airmass
    cmd = cmd + ' -s "' + offset
    return cmd

# Write an observing command prompting for filters, exptime, n_exp
# lunardist, airmass limits are set by default
# initial offset is always 1m and no dithering, which is not ideal
def make_obs_script_prompt(queue_id_num):
    # Don't use filter names yet, until order is verified.
    # filter_names = set_filter_dict()
    filter_names = None
    cmd = set_target_defaults()
    print 'Filter number is order in filter wheel GUI starting at 0; poss U,R,B,V,I,schott'
    exp_string = ''
    while True:
        obs_str = raw_input('Enter filter number, exptime, num. of exp [hit return to stop]: ')
        if obs_str.strip() == '':
            break
        fields = obs_str.split()
        if filter_names is not None:
            try:
                filtnum = filter_names[fields[0]]
            except:
                filtnum = int(fields[0])
        else:
            filtnum = int(fields[0])
        exptime = float(fields[1])
        nexp = int(fields[2])
        # Truncating exptime to an integer for now
        one_exp = ' filter=%1i E %i' % (filtnum, exptime)
        for i in range(nexp):
            exp_string = exp_string + one_exp
    # end loop
    cmd = cmd + exp_string
    queue_id_str = '%i' % (queue_id_num)
    cmd = cmd + '" ' + queue_id_str
    print "Command with script is:"
    print cmd
    return cmd

# Make the obs script given the queue id and a list of
# [filtnum1, exptime1, n_exp1, ...]
def make_obs_script(queue_id_num, explist, filter_names=None):
    cmd = set_target_defaults()
    exp_string = ''
    # Need 3 entries for each filter
    nfilts = int(len(explist) / 3)
    for i in range(nfilts):
        j = 0 + i*3
        if filter_names is not None:
            # here explist[j] is a string with the filter name, eg 'R' or '1'
            # use try/except to catch an undefined filter. The except will
            # fall back to a filter number, so if not an integer, it will still fail.
            try:
                filtnum = filter_names[explist[j]]
            except:
                filtnum = int(explist[j])
        else:
            filtnum = int(explist[j])
        exptime = float(explist[j+1])
        nexp = int(explist[j+2])
        # Truncating exptime to an integer for now
        one_exp = ' filter=%1i E %i' % (filtnum, exptime)
        for i in range(nexp):
            exp_string = exp_string + one_exp
    # end loop
    cmd = cmd + exp_string
    queue_id_str = '%i' % (queue_id_num)
    cmd = cmd + '" ' + queue_id_str
    return cmd        
    
# Add object number N to a queue
def add_object_queue(queue_name, queue_id_num):
    cmd = 'rts2-queue --queue %s %i' % (queue_name, queue_id_num)
    print cmd
    return cmd

# Run all commands: make the target, make its script, add to queue
# write the commands to a file
def make_obs_interactive(fname):
    queue_id_num = make_newtarget_prompt(fname)
    cmd = make_obs_script_prompt(queue_id_num)
    f = open(fname,'a')
    f.write(cmd+'\n')
    queue_name = 'plan'
    cmd2 = add_object_queue(queue_name, queue_id_num)
    f.write(cmd2+'\n')
    f.write('\n')
    f.close()
    print "File ",fname," now has rts2 commands appended"
    return

# Read a line of
#   name ra dec queue-id filter_num1 exptime1 n_exp1 filter_num2 exptime2 n_exp2 ...
# and write commands
def make_commands_fromfile(cline,filter_names=None):
    fields = cline.split()
    if len(fields) < 7:
        print "Line for ",fields[0]," doesn't have enough data"
        return 0,'',''
    name = fields[0]
    ra = fields[1]
    dec = fields[2]
    queue_id_num = int(fields[3])
    expstuff = fields[4:]
    cmd1 = make_newtarget_args(name,ra,dec,queue_id_num)
    cmd2 = make_obs_script(queue_id_num,expstuff,filter_names=filter_names)
    queue_name = 'plan'
    cmd3 = add_object_queue(queue_name, queue_id_num)
    return queue_id_num, cmd1, cmd2, cmd3

# Read the targets and info from a file
def make_obs_fromfile(inputname,fname):
    filters = set_filter_dict()
    infile = open(inputname,'r')
    f = open(fname,'a')
    for line in infile:
        # skip empty lines and comment lines
        if (line.strip() != '' and line[0] != '#'):
            # Don't use the filter dict until we confirm filter order.
            # queue_id_num, cmd1, cmd2, cmd3 = make_commands_fromfile(line,filter_names=filters)
            queue_id_num, cmd1, cmd2, cmd3 = make_commands_fromfile(line)
            f.write(cmd1+'\n')
            f.write(cmd2+'\n')
            # queue_name = 'plan'
            f.write(cmd3+'\n')
            f.write('\n')
    # end for loop
    infile.close()
    f.close()
    print "File ",fname," now has rts2 commands appended"
    return


def main():
    inputname = raw_input('Enter file with targets, or hit return to enter interactively: ')
    fname = raw_input('Enter filename to write rts2 commands: ')
    if inputname.strip() == '':
        make_obs_interactive(fname)
    else:
        print 'Target file format is: name ra dec queue-id filter_num1 exptime1 n_exp1 filter_num2 exptime2 n_exp2 ...'
        make_obs_fromfile(inputname,fname)
        
    print 'Execute the file with: sh ',fname,' , or similar'
    return

# This is the standard boilerplate that calls the main() function.
#if __name__ == '__main__':
  #main()

    
    
