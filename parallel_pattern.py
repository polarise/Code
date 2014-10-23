import multiprocessing

class Worker( multiprocessing.Process ):
	def __init__( self, i, q1, q2 ):
		"""
		i		- process_id
		q1	- input queue
		q2	- output queue
		"""
		multiprocessing.Process.__init__( self ) # inherit from ~
		self.i = i
		self.q1 = q1
		self.q2 = q2
	
	def fxn( self, myobject ):
		"""
		what each process should do
		"""
		print myobject.strip( "\n" ) # trivial operation
	
	def run( self ):
		"""
		overwrite the run method of the multiprocessing.Process class
		"""
		while True:
			myobject = self.q1.get()
			if myobject == None: # the poison pill; assume that 'None' is not a valid myobject value
				self.q2.put( "END" ) # mark that this process is done
				self.q1.task_done() # unblock q1
				self.q2.task_done() # unblock q2
				break # die
			else:
				result = self.fxn( myobject )
				self.q2.put( result )
				self.q1.task_done() # unblock q1
				self.q2.task_done() # unblock q2

if __name__ == "__main__":
	try:
		fn = sys.argv[1]
	except IndexError:
		print >> sys.stderr, "usage:./script.py <filename>"
		sys.exit( 1 )
	
	num_of_processors = multiprocession.cpu_count()
	q1 = multiprocessing.JoinableQueue() # input queue
	q2 = multiprocessing.JoinableQueue() # output queue
	
	# load the container e.g. from a file
	with open( fn ) as input_container:
		for myobject in input_container:
			q1.put( myobject )
	
	# put the poison pill
	for i in xrange( num_of_processors ):
		q1.put( None ) # 'None' is the poison pill
	
	# define the child processes, one process per processor
	procs = [ Worker( i, q1, q2 ) for i in xrange( num_of_processors )]
	
	# start the processes
	for p in procs:
		p.daemon = True # start as a daemon
		p.start()
		
	# block q1
	q1.join()
	
	# getting results
	end = list() # keep track of dead processes
	output_container = list() # get the results
	while len( end ) < num_of_processors:
		result = q2.get()
		if result != "END": # the end-of-processing signal
			output_container.append( result )
		else:
			end.append( result )
	
	# block q2
	q2.join()