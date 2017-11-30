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
#  BJW, Nov 29 2017
#
#  TODO: read input from file rather than interactive
#  TODO: have it execute the script rather than asking the user to run w/sh
#  TODO: have rts2-newtarget use command line args
#  TODO: dithering?

import os
import sys
# from astropy.io import ascii

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

# Write an observing command prompting for filters, exptime, n_exp
# lunardist, airmass limits are set by default
# initial offset is always 1m and no dithering, which is not ideal
def make_obs_script_prompt(queue_id_num):
    lunardist = '20:'
    airmass = ':2.2'
    camera = 'C0'
    offset = 'BIG61.OFFS=(1m,0)'
    cmd = 'rts2-target'
    cmd = cmd + ' -c ' + camera
    cmd = cmd + ' --lunarDistance ' + lunardist
    cmd = cmd + ' --airmass ' + airmass
    cmd = cmd + ' -s "' + offset
    print 'Filter number is order in filter wheel GUI starting at 0; poss U,R,B,V,I,schott'
    exp_string = ''
    while True:
        obs_str = raw_input('Enter filter number, exptime, num. of exp [hit return to stop]: ')
        if obs_str.strip() == '':
            break
        fields = obs_str.split()
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

def main():
    fname = raw_input('Enter filename to write rts2 commands: ')
    make_obs_interactive(fname)
    print 'Execute the file with: sh ',fname,' , or similar'
    return

# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
  main()

    
    
