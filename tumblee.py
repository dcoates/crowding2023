from psychopy import *
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


#Creates a GUI where the examiner enters the information of the subject that is to be tested
myDlg = gui.Dlg(title="Reading Experiment")
myDlg.addText('Subject Information')
myDlg.addField('Name:','TEST')
myDlg.addField('IGNORE',1)
myDlg.addField('Use Pupil Labs Eye Tracker', choices=["False", "True"])
myDlg.addField('Ecc',choices=["10", "20"])
myDlg.addField('Size (deg):',1)
myDlg.addField('Which letter:',choices=['E','T'])
myDlg.addField('Contrast:',choices=[1.0,0.1,0.01])
myDlg.addField('Duration',choices=[0.05,0.075,0.15,0.3,-1])
ok_data = myDlg.show()  # show dialog and wait for OK or Cancel
if myDlg.OK:  # or if ok_data is not None
    SubjectName = ok_data[0]
    subj_id=ok_data[1]
    use_tracker = ok_data[2]
    ecc = float(ok_data[3])
    size_deg = float(ok_data[4])
    which_let = ok_data[5]
    contrast = float(ok_data[6])
    trial_time = float(ok_data[7])
else:
    print('user cancelled')

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
targcol = [ 0.0 - contrast ] *3

fix_pos = ( (-1920/2*7/8,0) ) # Fixation in LVF 7/8 to the end of the screen

distance = 200
height_mm = 100

exper = conditions.experSony( (1920,1080), distance, height_mm, 9.0, 370 )
font = conditions.fontArial( targcol )
font.setCharDegs( stimHeightDeg, exper )

# This experiment - trials, timing, etc. 
ntrials = 1000.0
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


fixation = visual.TextStim(myWin,pos=(0,exper.yloc_fixation_pix),alignHoriz='center',height=9, color=(-1,-1,-1), ori=0, font=font.selfont,
                           anchorHoriz='center',anchorVert='center' ) 
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

size_pix = exper.deg2pix( size_deg ) # +font.let_height_ptfont/2.0
targ = stims.stim_letter( myWin, size_pix, font.contrast, font.selfont, 
	{'height': test_heights, 'xpos':0.0, 'ypos':0, 'text': which_let, 'ori':targseq}  )
left = stims.stim_letter( myWin, size_pix, font.contrast, font.selfont, 
	{'xpos':-spacing, 'ypos':exper.yloc_pix, 'text': which_let, 'height':test_heights, 'ori':seqL,'anchorHoriz':'center', 'anchorVert':'center'}  )
right = stims.stim_letter( myWin, size_pix, font.contrast, font.selfont, 
	{'xpos':+spacing, 'ypos':exper.yloc_pix, 'text': which_let, 'height':test_heights, 'ori':seqR,'anchorHoriz':'center', 'anchorVert':'center'}  )
up = stims.stim_letter( myWin, size_pix, font.contrast, font.selfont, 
	{'xpos':+spacing, 'ypos':exper.yloc_pix, 'text': which_let, 'height':test_heights, 'ori':seqU,'anchorHoriz':'center', 'anchorVert':'center'}  )
down = stims.stim_letter( myWin, size_pix, font.contrast, font.selfont, 
	{'xpos':+spacing, 'ypos':exper.yloc_pix, 'text': which_let, 'height':test_heights, 'ori':seqD,'anchorHoriz':'center', 'anchorVert':'center'}  )

stims = [ targ, right, up, left, down, ] #, olR]


for stim1 in stims:
    stim1.widget.anchorHoriz='center'
    stim1.widget.anchorVert='center'

for stim1 in [fixation]:
    stim1.anchorHoriz='center'
    stim1.anchorVert='center'

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
msg = 'Calibrating monitor timing...'
fixation.setText( msg.upper() ) # sloan needs to be uppercase
for i in numpy.arange(20):
	fixation.draw()
	[stim.draw() for stim in stims]
	#[cue.draw() for cue in cues]
	myWin.flip()
savetimes = myWin.frameIntervals
fliprate = numpy.mean( savetimes[-40:] )
print( 'fliprate=%f' % fliprate)

fixation.draw()
myWin.flip()
event.waitKeys()

fixation.setPos( fix_pos) # Fixation in LVF 7/8 to the end of the screen
fixation.setHeight(200)
fixation.setText( '+' )
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

	if repeat == False:
		[stim.getTrial(trialNum) for stim in stims]

	try:
		sval = thisStair.next()
	except StopIteration:#we got a StopIteration error
		print ('Finished. intensity: %f' % (numpy.mean( thisStair.reversalIntensities[2:] ) ) )
		break

#	print sval
	#font.setCharDegs( sval, exper )
	#for stim in stims:
		#stim.height=font.let_height_ptfont

	cx = fix_pos[0] + exper.deg2pix(ecc)
	cy = 0
	targ.pos =  (cx,cy)
	left.pos =  (cx - sval * size_pix, cy )
	right.pos = (cx + sval * size_pix, cy )
	up.pos =    (cx, cy + sval * size_pix )
	down.pos =  (cx, cy - sval * size_pix )

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

	valid_resp = True

	if which_let=='E':
		resp_map={'left': 180, 'right':0, 'up':270, 'down':90}
	elif which_let=='T':
		resp_map={'left': 270, 'right':90, 'up':0, 'down':180}
	for key in event.waitKeys():
		if key in [ 'escape', 'q' ]:
			done = True
			valid_resp = False
		else:
			resp_ori = resp_map[key]

	iscorrect = (resp_ori == targ.ori)
	print(key, targ.ori, iscorrect, valid_resp,resp_ori,targ.strvals())
	if valid_resp:
	    iscorrect = (resp_ori == targ.ori)
	    s="%f,%s,"%(sval,key)
	    s+= ",".join([str(s.ori) for s in stims])
	    s+= ',%d'%int(iscorrect)
	    outfile.writelines('%s\n'%s )

	    thisStair.addData( iscorrect )

	    trialNum += 1
	    if trialNum >= maxtrials:
		    done = True

myWin.close()

outfile.write( 'height=' + str(stimHeightDeg) + "\n")
outfile.write ('mean value: %f\t' % (numpy.mean( thisStair.reversalIntensities[2:] ) ) )
outfile.write ('SubjectName: %s\n' %  SubjectName )

if OutputHeader:
	outfile.write( '#END\n' )
	outfile.write( "#V 0.2\n" )
	outfile.write( '#contrast=' + str(contrast) + "\n" )
	outfile.write( '#size_pix=' + str(size_pix) + "\n")
	outfile.write( '#size_deg=' + str(size_deg) + "\n")
	outfile.write( '#ecc=' + str(ecc) + "\n")
	outfile.write( '#distance=' + str(distance) + "\n")
	outfile.write( '#cx=' + str(cx) + "\n" )
	outfile.write( '#cy=' + str(cy) + "\n" )
	outfile.write( '#trial time=' + str(trial_time) + "\n")
	outfile.write( '#which_let' + str(which_let) + "\n")

outfile.close()
