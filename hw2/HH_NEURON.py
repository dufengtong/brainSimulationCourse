# Copyright 2018 Du Fengtong
#
# PyNeuron implication of HH neuron with different current clamp.
# four kinds of current clamp:
# 1. direct current
# 2. direct current with gaussian noise
# 3. sinusoidal current
# 4. ramp current

from neuron import h
from pylab import *

# run control
h.load_file("stdrun.hoc")

# set membrane property
cell = h.Section()
cell.insert('hh')

# record voltage and time
vrec = h.Vector()
trec = h.Vector()
vrec.record(cell(0.5)._ref_v) # 0.5 means the position of the recorder
trec.record(h._ref_t)

# set the running duration and time step
h.tstop = 1e3
h.dt = 0.025

# set current clamp type
clamp_type = 'noise'

# set current clamp position and duration
icc = h.IClamp(cell(0.5))
icc.delay=0
icc.dur=1e9

if clamp_type == 'noise':
    AMP = 20/1.5
    np.random.seed(21051982)
    noise = AMP*(1+np.random.rand(int(h.tstop/h.dt+1)))
    nvec = h.Vector( noise )
    nvec.play(icc._ref_amp,h.dt)
elif clamp_type == 'direct':
    AMP = 20
    inj = np.ones(int(h.tstop/h.dt+1))*AMP
    nvec = h.Vector( inj )
    nvec.play(icc._ref_amp,h.dt)
elif clamp_type == 'ramp':
    AMP = 25
    inj = 16 + np.linspace(0,AMP,num=int(h.tstop/h.dt+1))
    nvec = h.Vector( inj )
    nvec.play(icc._ref_amp,h.dt)
elif clamp_type == 'sin':
    AMP = 20
    t = np.linspace(0, h.tstop, num=int(h.tstop/h.dt+1))
    inj = AMP*(sin(0.1*t)+1)
    nvec = h.Vector( inj )
    nvec.play(icc._ref_amp,h.dt)
else:
    raise Exception("Unsupported clamp type, please choose from "
                    "'noise', 'direct', 'ramp' and 'sin'")

# initial membrane voltage
h.v_init = -60
# start running
h.run()

# plot results
subplot(2,1,1)
a = nvec.as_numpy()
plot(trec.as_numpy(),nvec.as_numpy()) # display part of the noisy current-injection amplitude
xlim((0,h.tstop))
xlabel('Time (ms)')
ylabel('i (nA)')
title('inject current')
ylim([0,40])

subplot(2,1,2)
plot(trec.as_numpy(),vrec.as_numpy()) # as_numpy is more efficient than converting/copying Vector to numpy format
xlim((0,h.tstop))
xlabel('Time (ms)')
ylabel('Vm (mV)')
title('HH neuron')
tight_layout()

show()
