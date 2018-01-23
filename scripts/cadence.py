from rts2 import scriptcomm
import datetime
import rts2


class Cadence(scriptcomm.Rts2Comm):
	"""First attempt at making a script for a cadenced observation"""

	def __init__(self, start_times, start_window, target_id, qname="plan" ):
		"""Constructor

        Args:
            start_times => a list of datetimes which are the start times for 
            a cadenced observation

            start_window => the longest amount of time the object can be usefully observed 
            after the start time. 


            Example:
                If we wanted to observe an object at 19:30, 21:30 and 23:30
                we could set the start times to be 19:25, 21:25 and 23:25
                respectively and set the start_window to 5 minutes. 
                
        """

		
		self.start_time = start_times

		#The longest time you would be willing to 
		#to let the observation start (after start_time)
		self.start_window = start_window

		self.prx = rts2.json.createProxy('http://localhost:8889',username='petr',password='test')
        self.target_id = target_id
        self.qname = qname


	def run( self ):
        now = datetime.datetime.now()
        name_format = "%b/%N.fits"
        for start_time in self.start_times:
            if (now > start_time and ( start_time + self.start_window ) ):
                #We are within the window expose!
                img = self.expose(self.before_readout, 10, name_format )
            
        else
            # we are not in the window stick this observation later in the queue. 
            q=rts2.queue.Queue(self.prx, "manual")
			q.load()
			
			#remove this target id from queue
			entries = [e.target.id for e in q.entries if e.target.id is not self.target_id]
			
			#put it back after the next observation
			entries.insert( 1, self.target_id )
			
				
			
			


    def before_readout(self):
        #set filters and stuff here
        pass



if __name__ == "__main__":
	now = datetime.datetime.now()
	obs1 = datetime.datetime(now.year, now.month, now.day, 19, 25, 00)
	obs2 = datetime.datetime(now.year, now.month, now.day, 21, 25, 00)
	obs3 = datetime.datetime(now.year, now.month, now.day, 23, 25, 00)
	window = datetime.timedelta(minutes=6)
	obs = Cadence( [obs1, obs2, obs3], window, 636 )
	obs.run()
