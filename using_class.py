class Gene:
	def __init__( self, name ):
		self.name = name
	
	def say_my_name( self ):
		print self.name
	
	def say_name_in_swahili( self, someone ):
		print "%s, jina langu ni %s." % ( someone, self.name )
	
G = Gene( "HIST34" )
G.say_my_name()
G.say_name_in_swahili( "Claire" )
