# Copyright 2018 Du Fengtong
#
# simulate regular spiking with adaptation

from neuron import h
from pylab import *

# run control
h.load_file("stdrun.hoc")

# set membrane property
cell = h.Section()
cell.insert('burste')

# record voltage and time
vrec = h.Vector()
trec = h.Vector()
vrec.record(cell(0.5)._ref_v) # 0.5 means the position of the recorder
trec.record(h._ref_t)

# set the running duration and time step
h.tstop = 1e3
h.dt = 0.025

# set current clamp position and duration
icc = h.IClamp(cell(0.5))
icc.delay=200
icc.dur=600
AMP = 0.7 # nA
inj = np.zeros(int(h.tstop/h.dt+1))
inj[int(icc.delay/h.dt):int((icc.delay+icc.dur)/h.dt)] = AMP
irec = h.Vector(inj)
irec.play(icc._ref_amp,h.dt)

# initial membrane voltage
# h.v_init = -60
# start running
h.run()

# plot results
subplot(2,1,1)
plot(trec.as_numpy(),irec.as_numpy()) # display part of the noisy current-injection amplitude
xlim((0,h.tstop))
xlabel('Time (ms)')
ylabel('i (nA)')
title('inject current')
# ylim([0,40])

subplot(2,1,2)
plot(trec.as_numpy(),vrec.as_numpy()) # as_numpy is more efficient than converting/copying Vector to numpy format
xlim((0,h.tstop))
xlabel('Time (ms)')
ylabel('Vm (mV)')
title('Bursting-inhibitory Neuron')
tight_layout()
show()
