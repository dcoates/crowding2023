import psychopy, numpy
import os
  
subj = ' yytyfilt20'
#subj = 'xx'
backs = [ 14 ] # make 16,3,20,then 1

#spacings = [1.5, 1.625, 1.75, 2.0, 4.0, 99.0 ]
spacings = [1.5, 3.0, 99.0 ]
#spacings = [1.5, 1.75, 4.0, 99.0 ]
#spacings = [ 1.625,  2.0 ]
#spacings = [ 1.75 ]
#spacings = [ 2.0 ]
#spacings = [  ]
#spacings = [ 99 ]
for aback in backs:
	allspacings = numpy.tile( spacings, (1, 2 ))[0]
	allspacings = numpy.random.permutation( allspacings )
	print allspacings
	
	exename = "c:\Progra~2\PsychoPy2\python.exe" # Windows
#	exename = "python" # Linux
	
	for aspacing in allspacings:
		fname = "%s-%s-%s" % (subj, aback, str(aspacing))
		print allspacings, aspacing
		os.system('%s tumblee.py %s %s' % (exename, fname,str(aspacing)) )
