import psychopy
import psychopy.visual as visual
import psychopy.event as event
import psychopy.core as core
import numpy as np
import time

if __name__ == "__main__":

    done=False

    win = visual.Window(allowGUI=True, units='pix', size=(1024,768), fullscr=True,color=[-1,-1,-1])
    txt=visual.TextStim(win,"CON", pos=(500,300), height=100.0, color=[0,1,1] )
    win.setMouseVisible(False)

    stims=[txt]

    con=0.0

    while done==False:

        if con<0: con=0
        if con>1: con=1
        txt.text='%0.1f'%con
        win.color = [con*2-1]*3

        [astim.draw() for astim in stims]
        win.flip()
        keys=event.waitKeys()

        for key in keys:
            if key in [ 'escape', 'q' ]:
                done = True
            if key in ['down','right']:
                con += 0.1
            elif key in ['up','left']:
                con -= 0.1
