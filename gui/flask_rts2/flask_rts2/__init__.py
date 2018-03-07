from __future__ import print_function
from flask import Flask, send_file
from flask import render_template
from flask_basicauth import BasicAuth
import requests
import json
import glob
import os
import sys
#from rts2 import rtsapi

app = Flask(__name__)
app.config.from_object(__name__)

CONFIG_PATH = "/home/rts2obs/.mtnops"
CONFIG_FILE = "flask_rts2.json"
CONFIG_FD = open( os.path.join( CONFIG_PATH, CONFIG_FILE ) )

CONFIG = json.load(CONFIG_FD)

print( CONFIG, file=sys.stderr ) 
CONFIG_FD.close()



app.config['BASIC_AUTH_USERNAME'] = CONFIG["username"]
app.config['BASIC_AUTH_PASSWORD'] = CONFIG['password']

basic_auth = BasicAuth( app )



# Begin worker fucntions

def _get_device( device ):
    """
    args: device -> name of the device

    Description: Uses RTS2 HTTP interface to get a json object of all the rts2 values 
        (the ones you see in rts2-mon) of a device
    
    returns json with device data or error message
    """
    try:
        r=requests.get( "http://localhost:8889/api/get?e=1&d={}".format( device ) )
        data = r.text
    except Exception as err:
        data = json.dumps( { "error": str(err) } )

    return data


def _set_rts2_value(device, name, value):
    """
    args: 
        device -> name of the device
        name -> name of the value to be set ie queue_only
        value -> value to set it to 

    Description:    uses rts2 http interface to set an RTS2 value.
    
    returns json with device data or error message
    """

    try:
        r=requests.get("http://localhost:8889/api/set?async=None&v={}&d={}&n={}".format(value, device, name))
        data=r.text
    except Exception as err:
        data = json.dumps({"error": err})

    return data


#end worker functions


#Begin flask functions

@app.route( '/' )
@basic_auth.required
def root():
    """Main page renders the index.html page"""
    return render_template( "index.html", username=app.config['BASIC_AUTH_USERNAME'], passwd=app.config['BASIC_AUTH_PASSWORD'] )



@app.route( '/device/<name>' )
@basic_auth.required
def get_device( name ):
    return _get_device(name)




@app.route('/device')
@basic_auth.required
def get_all_devices():
    """Get the rts2 values for all the devices in RTS2 and put them in
    one big json."""
    outdata = {}
    for dev in ( "BIG61", "C0", "SEL", "EXEC", "F0", "W0" ):
        jdata = _get_device( dev )
        try:
            jdata = json.loads( jdata )
        except Exception as err:
            jdata = {"error":str(err)}

        outdata[dev] = jdata

    return json.dumps(outdata)



@app.route('/device/set/<device>/<name>/<value>')
@basic_auth.required
def set_rts2_value(device, name, value):
    _set_rts2_value(device, name, value)



@app.route('/queuestart')
@basic_auth.required
def rts2_queue_start():
    """Easy way set all the parameters for starting
    queue observing. """
    _set_rts2_value("SEL", "plan_queing", 3)
    _set_rts2_value("SEL", "queue_only", True)
    _set_rts2_value("BIG61", "pec_state", 1)
    _set_rts2_value("EXEC", "auto_loop", False)
    _set_rts2_value("BIG61", "dome_auto", True)


@app.route('/lastimg')
@basic_auth.required
def download_lastimg():
    """Use C0.last_img_path to downlaod the most recent image."""
    jdata = _get_device("C0")
    return send_file( json.loads(jdata.text)["d"]["last_img_path"][1] )



@app.route('/weather/boltwood.json')
@basic_auth.required
def boltwood_json():
    """Use C0.last_img_path to downlaod the most recent image."""
    try:
        r=requests.get("https://www.lpl.arizona.edu/~css/bigelow/boltwoodlast.json")
        jdata = r.text
    except Exception as err:
        jdata = json.dumps({"error": str(err)})
    return jdata


if __name__ == "__main__":
    app.run(host="0.0.0.0")



