from psychopy import *
#from psychopy import log
import psychopy.data as pd 
import numpy
import time
import matplotlib.pyplot as pyplot
import os 
import pylab # for frame interval plotting only
import sys

import stims
import conditions

if len(sys.argv) == 3:
    SubjectName = sys.argv[1]     
    spacingMult = float(sys.argv[2])
else:
    SubjectName = 'xx'
    spacingMult = 4.0
    
stimHeightDeg = 2.25 # dummy

# Utility functions
def randi( high, numvals=None, low=0 ):
	if numvals == None:
		return numpy.array(numpy.random.rand() * (high-low)+low, dtype=int)
	else:
		return numpy.array(numpy.random.rand(numvals) * (high-low)+low, dtype=int)

def buildtrialseq( vals, ntrials ):
	return numpy.random.permutation(
        numpy.tile ( vals, int(numpy.ceil( float(ntrials)/len(vals)) ))
    )[0:int(ntrials)]

#s Outputfile params
CreateUniqueBlockFiles = True
OutputHeader = True

# This experiment: contrasts
#backcol= ( -0.5, -0.5, -0.5 )
backcol= ( 0,0,0 )
#targcol= (-1.00,-1.00,-1.00)	
targcol= ( -1.0, -1.0,  -1.0 )
#targcol= ( 1.0, 1.0, 1.0 )

# TODO: We want it at 8 degress, so why need 9??
exper = conditions.experSony( (1920,1080), 1080.0, 9.0, 370 )
font = conditions.fontArial( targcol )
font.setCharDegs( stimHeightDeg, exper )

# This experiment - trials, timing, etc. 
ntrials = 1000.0
trial_time =  0.150  # in sec. use -1 for infinite
pre_time = 0.00
mask_time = 0.000

trial_cycles = 0.0

fullscr = True # make sure this matches next

# Set up the screen, etc.
myWin = visual.Window(exper.screendim, allowGUI=True, color=backcol, units='pix', fullscr=fullscr) #, winType='pygame' )
myWin.setMouseVisible(False)
myWin.setRecordFrameIntervals(True)
#myWin._refreshThreshold=0.03 # set to 30 ms on my monitor
#mywin._refreshThreshold=1/85.0+0.004 #i've got 85Hz monitor and want to allow 4ms tolerance

#log.console.setLevel(log.ERROR)

fixation = visual.TextStim(myWin,pos=(0,exper.yloc_fixation_pix),alignHoriz='center',height=9, color=font.contrast, ori=0, font=font.selfont ) 
fixation.setText( 'o' )

# Post-trial text
response_list_disp = visual.TextStim(myWin,pos=(0,exper.yloc_fixation_pix),alignHoriz='center',height=font.let_height_ptfont,
	color=font.contrast, ori=0, font=font.selfont, text=' ' ) 

orda = 65 #97 
#targets = ['b', 'o', 'l', 'd', 'c', 'p', 'q', 'h', 'k', 'x', 'z', 's', 'e', 'a', 'u', 'f', 'g' ]
#targets = ['b', 'o', 'l', 'd', 'p', 'q' ]
#targets = ['b', 'd'  ]
#targets = ['b', 'd', 'o', 'l'  ]
#targets = [ 'b', 'd', 'o', 'l' ]
#targets = [ chr(i+orda) for i in numpy.arange(26) ]
targets = [0,90,180,270]
#targets = [ 'b', 'd', 'o', 'l' ]

#response_list_disp.setText( str(targets) )
response_list_disp.setText( "" )

maxtrials = ntrials
#targseq = [targets[i] for i in randi(len(targets), maxtrials) ]
test_heights = font.let_height_ptfont

spacseq = buildtrialseq( [exper.deg2pix(angl) for angl in [1.5]], ntrials )
spacing = exper.deg2pix( stimHeightDeg ) # +font.let_height_ptfont/2.0

targseq = [targets[i] for i in buildtrialseq( numpy.arange(len(targets)), maxtrials) ]
seqL = [targets[i] for i in buildtrialseq( numpy.arange(len(targets)), maxtrials) ]
seqR = [targets[i] for i in buildtrialseq( numpy.arange(len(targets)), maxtrials) ]
seqU = [targets[i] for i in buildtrialseq( numpy.arange(len(targets)), maxtrials) ]
seqD = [targets[i] for i in buildtrialseq( numpy.arange(len(targets)), maxtrials) ]
 
targ = stims.stim_letter( myWin, font.let_height_ptfont, font.contrast, font.selfont, 
	{'height': test_heights, 'xpos':0.0, 'ypos':exper.yloc_pix, 'text':'E', 'ori':targseq}  )
left = stims.stim_letter( myWin, font.let_height_ptfont, font.contrast, font.selfont, 
	{'xpos':-spacing, 'ypos':exper.yloc_pix, 'text':'E', 'height':test_heights, 'ori':seqL}  )
right = stims.stim_letter( myWin, font.let_height_ptfont, font.contrast, font.selfont, 
	{'xpos':+spacing, 'ypos':exper.yloc_pix, 'text':'E', 'height':test_heights, 'ori':seqR}  )

up = stims.stim_letter( myWin, font.let_height_ptfont, font.contrast, font.selfont, 
	{'xpos':+spacing, 'ypos':exper.yloc_pix, 'text':'E', 'height':test_heights, 'ori':seqU}  )
down = stims.stim_letter( myWin, font.let_height_ptfont, font.contrast, font.selfont, 
	{'xpos':+spacing, 'ypos':exper.yloc_pix, 'text':'E', 'height':test_heights, 'ori':seqD}  )

if spacingMult < 99:
	stims = [ left, targ, right  ] #, olR]
else:
	stims = [ targ  ] #, olR]

#stims = [ otarg ]
#stims = [ ]

trial_cyles = 1

cues = [] # [ cueU, cueD ]

# Noise mask
noise_size = 128  # TODO: Parametrize. Base on min/mix of stimuli??
noise = visual.PatchStim( myWin, texRes=1, mask="none", tex="none", pos=(0,-exper.deg2pix(7.5)), units='pix', size=(noise_size, noise_size), color=[1.0,1.0, 1.0] )
noisemult = 2.0 # multiplier for rand (contrast)
noiseoff =  2.0 # - offset of rand
noisebinary = True # TODO: Better specify contrast of binary noise: is floor(rand*mult)-off

maxtrials = ntrials

if CreateUniqueBlockFiles:
	# Find a unique filename by appending a number to the session
	blockidx = 0
	gotunique = False
	while gotunique==False: 
		outfilename = "%s_%s-%02d.csv" % (SubjectName, time.strftime("%m%d%Y", time.localtime() ), blockidx)
		if os.path.exists(outfilename):
			blockidx += 1
		else:
			gotunique = True
else:
	outfilename = "%s_%s.csv" % (SubjectName, time.strftime("%m%d%Y", time.localtime() ) )
		
outfile = open(outfilename, "wt")
xval = 0
done = False
trialNum=0

# Calibrate by seeing how long 100 redraws takes
fixation.setHeight(30)
fixation.setPos( (0,0) )
msg = 'Hello! %s -- Spacing is: %f -- Please note filter-- Thx!!' %( SubjectName, spacingMult )
fixation.setText( msg.upper() ) # sloan needs to be uppercase
#fixation.setText( 'Calibrating monitor...' )
for i in numpy.arange(100):
	fixation.draw()
	[stim.draw() for stim in stims]
	#[cue.draw() for cue in cues]
	myWin.flip()
savetimes = myWin.frameIntervals
fliprate = numpy.mean( savetimes[20:80] )
print( 'fliprate=%f' % fliprate)

#ixation.setHeight(20)
#sg = 'Hello! %s\n Spacing is: %f\n Please remember filter/projector.\nThanks!' %( SubjectName, spacingMult )
#ixation.setText( msg.upper() ) # sloan needs to be uppercase
#fixation.setText( 'Press any key twice to start (first goes to fixation screen)\n'.upper() )
fixation.draw()
myWin.flip()
event.waitKeys()
fixation.setPos( (0,exper.yloc_fixation_pix) )

fixation.setHeight(70)
fixation.setText( '-----' )
fixation.draw()
myWin.flip()
event.waitKeys()

repeat = False

# TODO: was 5
thisStair = pd.StairHandler(startVal=4, nTrials=50, nUp=1, nDown=3, minVal = 0.5, maxVal=7, stepSizes=[3,2,1,0.75,0.5,0.5,0.2,0.2,0.2,0.2,0.2,0.2] ) #, stepSizes=[4,2,1,1,1,1,1,1])
#thisStair = pd.StairHandler(startVal=50, nTrials=50, nUp=1, nDown=3, minVal = 0.5, maxVal=50, stepSizes=[3,2,1,0.75,0.5,0.5,0.2,0.2,0.2,0.2,0.2,0.2] ) #, stepSizes=[4,2,1,1,1,1,1,1])

#thisStair = pd.StairHandler(startVal=50, nTrials=50, nUp=1, nDown=3, minVal = 0.0, maxVal=100, stepType='lin', stepSizes=[10.,5,1,1,1,1,1,1,1] ) #, stepSizes=[4,2,1,1,1,1,1,1])

while not done:
	fixation.draw()
	myWin.flip()
	core.wait(pre_time)
	fixation.draw()

	#if (numpy.random.rand() < 0.5) & (trialNum>0):
		#repeat = True
	#else:
		#repeat = False
	
	if repeat == False:
		[stim.getTrial(trialNum) for stim in stims]

	try:
		sval = thisStair.next()
	except StopIteration:#we got a StopIteration error
		print ('Finished. intensity: %f' % (numpy.mean( thisStair.reversalIntensities[2:] ) ) )
		break

#	print sval
	font.setCharDegs( sval, exper )
	for stim in stims:
		stim.height=font.let_height_ptfont
#	left.pos =  (spacingMult * exper.deg2pix(sval), exper.yloc_pix )
#	right.pos =  (-spacingMult * exper.deg2pix(sval), exper.yloc_pix )
	left.pos =  (spacingMult * exper.deg2pix(sval), exper.yloc_pix )
	right.pos =  (-spacingMult * exper.deg2pix(sval), exper.yloc_pix )
#	left.pos =  ( exper.deg2pix ( spacingMult ), exper.yloc_pix )
#	right.pos =  (- exper.deg2pix(spacingMult), exper.yloc_pix )

	myWin.flip()
#	targcol = ( (sval-100)/100.0, (sval-100)/100.0, (sval-100)/100.0 )
	if trial_cycles > 0:
		numflips = trial_time / fliprate
		for angle in numpy.linspace( 0, numpy.pi*(trial_cycles/1.0), numflips ):
			#contr = -0.5 + -0.5 * abs(numpy.sin( angle) )
           
			[stim.draw( (contr,contr,contr) ) for stim in stims]
			[cue.draw() for cue in cues]
			fixation.draw()
			myWin.flip()
	else:
			[stim.draw( targcol ) for stim in stims]
			[cue.draw() for cue in cues]
			fixation.draw()
			myWin.flip()
			if ( trial_time >= 0  ):
				core.wait(trial_time)

	trialNum += 1
	#[stim.getTrial(trialNum) for stim in stims]
	if trial_cycles > 0:
		targ.getTrial(trialNum)
		numflips = trial_time / fliprate
		for angle in numpy.linspace( 0, numpy.pi*(trial_cycles/1.0), numflips ):
			contr = -0.5 + -0.5 * abs(numpy.sin( angle) )
			[stim.draw( (contr,contr,contr) ) for stim in stims]
			[cue.draw() for cue in cues]
			fixation.draw()
			myWin.flip()

	if mask_time > 0:
		fixation.draw()
		#response_list_disp.draw()

#		pct_value = randi(4)+1
#		noisemult = (1.0 + 1.0/(1.0/(0.1*pct_value)-1.0))
#
#		if numpy.random.rand() < 0.5:
#			is50 = True
#			pct_value += 5
#		else:
#			is50 = False
		
		if noisebinary:
			thenoise = numpy.floor(numpy.random.rand(noise_size,noise_size)*noisemult)-noiseoff 
#			if is50:
#				thenoise = (thenoise==0)*-1.0
			noise.setTex( thenoise )
		else:
			noise.setTex( (numpy.random.rand(noise_size,noise_size)*noisemult)-noiseoff)

		noise.draw()
		myWin.flip()
		core.wait(mask_time)

	# Allow infinite display
	if (trial_time >= 0):
		fixation.draw()
		response_list_disp.draw()
		myWin.flip()

	resp_ori = numpy.floor(numpy.random.rand()*4)*90
	for key in event.waitKeys():
		if key in [ 'escape' ]:
			#core.quit()
			done = True
		if key in [ 'left' ]:
			resp_ori = 180
		if key in [ 'right' ]:
			resp_ori = 0
		if key in [ 'up' ]:
			resp_ori = 90
		if key in [ 'down' ]:
			resp_ori = 270 
		if key in [ 'r' ]:
			repeat = (repeat == False)

	#outfile.write( "%s, %s %s %s\n" % (key, targ.strvals(), olL.strvals(), olR.strvals() ) )
	iscorrect = (resp_ori == targ.ori)
	thisStair.addData( iscorrect )

	trialNum += 1
	if trialNum >= maxtrials:
		done = True

myWin.close()

outfile.write( 'height=' + str(stimHeightDeg) + "\n")
outfile.write ('mean value: %f\t' % (numpy.mean( thisStair.reversalIntensities[2:] ) ) )
outfile.write ('spacingMult: %s\t' % str( spacingMult ) )
outfile.write ('SubjectName: %s\n' %  SubjectName )

if OutputHeader:
	outfile.write( '#END\n' )
	outfile.write( '%s\n' % [str(stim) for stim in stims])
	outfile.write( "#V 0.2\n" )
	outfile.write( '#font=' + str(font.let_height_ptfont) + "\n")
	outfile.write( '#font.let_height_pixels=' + str(font.let_height_pixels) + "\n")
	outfile.write('#font.let_height_mm=' + str(font.let_height_mm) + "\n")
	outfile.write('#font.let_height_deg=' + str(font.let_height_deg) + "\n")
	#outfile.write('#ratio of (ascender+descender) to o=' + str(ascender_and_descender_to_o) + "\n")
	#outfile.write('#ratio of ptfont to pixels=0' + str(ptfont_to_pixels) + "\n")
	#utfile.write( '#distance=' + str(distance) + "\n")
	#utfile.write( '#pixels per mm=' + str(screensize) + "\n")
	#utfile.write( '#degrees_eccentricity=' + str(degrees_eccentricity) + "\n" )
	outfile.write( '#exper.yloc_pix=' + str(exper.yloc_pix) + "\n" )
	outfile.write( '#trial time=' + str(trial_time) + "\n")
	outfile.write( "key, targ.strvals(), olL.strvals(), olR.strvals() ")

outfile.close()

#pylab.plot(myWin.frameIntervals)
#pylab.grid()
#pylab.show()


