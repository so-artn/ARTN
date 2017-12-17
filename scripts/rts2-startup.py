#!/usr/bin/python


"""This is the RTS2 startup script
this should be run at the beginning of the 
night.
    
    """

import sys
import rts2.rtsapi


if len(sys.argv) != 3:
    print "usage: {} <username> <password>".format(sys.argv[0])
    sys.exit()
username, password = sys.argv[1:]

prxy = rts2.rtsapi.createProxy("http://bigpop:8889", username=username, password=password)


prxy.setValue( "SEL", "queue_only", 1 )

#sets to WESTEAST
prxy.setValue( "SEL", "plan_queueing", 3 )
prxy.setValue( "BIG61", "pec_state", 1 )
prxy.setvalue( "EXEC", "auto_loop", False )


