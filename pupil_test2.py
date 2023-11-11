import psychopy
import psychopy.visual as visual
import psychopy.event as event
import psychopy.core as core
import numpy as np
import time
import sys

c_g=[0,1,0]
c_k=[-1,-1,-1]

# Use defaults or read reps & duration from cmd line
reps=2
secs_per=5
if len(sys.argv)>2:
    reps=int(sys.argv[1])
    secs_per=int(sys.argv[2])

# See if we are on a computer with the Pupil Labs API
if 1:
    from pupil_labs.realtime_api.simple import discover_one_device
    device = discover_one_device()
    USE_PUPIL=True
    print(f"Phone IP address: {device.phone_ip}")
    print(f"Phone name: {device.phone_name}")
    print(f"Battery level: {device.battery_level_percent}%")
    print(f"Free storage: {device.memory_num_free_bytes / 1024**3:.1f} GB")
    print(f"Serial number of connected glasses: {device.module_serial}")

    recording_id = device.recording_start()
    print(f"Started recording with id {recording_id}")
    device.send_event("START")
else:
    USE_PUPIL=False

if __name__ == "__main__":

    win = visual.Window(allowGUI=True, units='pix', size=(1024,768), fullscr=False,color=[0,0,0])
    win.setMouseVisible(False)

    fix=visual.TextStim(win,pos=(-300,0),anchorHoriz='center', anchorVert='center', text="+",color=[-1,-1,-1],height=100 )
    target=visual.TextStim(win,pos=(300,0),anchorHoriz='center', anchorVert='center', text="+",color=[-1,-1,-1],height=100 )

    message=visual.TextStim(win,pos=(0,-300),anchorHoriz='center', anchorVert='center', text="---",color=[-1,-1,-1],height=32 )
    if USE_PUPIL:
        message.text = "Tracking OK"
    else:
        message.text = "No tracking"

    fix.autoDraw = True
    target.autoDraw = True
    message.autoDraw=True

    nrep=0
    while nrep<reps:
        fix.color=c_g
        target.color=c_k
        if USE_PUPIL:
            device.send_event("%d_left"%nrep)
        win.flip()

        core.wait(secs_per)

        fix.color=c_k
        target.color=c_g
        if USE_PUPIL:
            device.send_event("%d_right"%nrep)
        win.flip()

        core.wait(secs_per)
        nrep += 1

    if USE_PUPIL:
        device.recording_stop_and_save()
