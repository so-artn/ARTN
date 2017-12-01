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
#  BJW, Nov 29-30 2017
#
#  DONE: read input from file rather than interactive
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
    cmd = set_target_defaults()
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

# Make the obs script given the queue id and a list of
# [filtnum1, exptime1, n_exp1, ...]
def make_obs_script(queue_id_num, explist):
    cmd = set_target_defaults()
    exp_string = ''
    # Need 3 entries for each filter
    nfilts = int(len(explist) / 3)
    for i in range(nfilts):
        j = 0 + i*3
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
def make_commands_fromfile(cline):
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
    cmd2 = make_obs_script(queue_id_num,expstuff)
    queue_name = 'plan'
    cmd3 = add_object_queue(queue_name, queue_id_num)
    return queue_id_num, cmd1, cmd2, cmd3

# Read the targets and info from a file
def make_obs_fromfile(inputname,fname):
    infile = open(inputname,'r')
    f = open(fname,'a')
    for line in infile:
        # skip empty lines and comment lines
        if (line.strip() != '' and line[0] != '#'):
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
if __name__ == '__main__':
  main()

    
    
