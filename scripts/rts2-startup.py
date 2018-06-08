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

#create the RTS2 proxy.
#This is a JSON http based interface
#for communication with rts2. 
prxy = rts2.rtsapi.createProxy("http://localhost:8889", username=username, password=password )

prxy.setValue( "SEL", "queue_only", True )

#sets to WESTEAST
prxy.setValue( "SEL", "plan_queing", 3 )
prxy.setValue( "BIG61", "pec_state", True )
prxy.setValue( "EXEC", "auto_loop", False )


