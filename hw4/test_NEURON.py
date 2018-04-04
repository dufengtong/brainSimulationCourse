from neuron import h
from pylab import *
import timeit

start = timeit.default_timer()


cell = h.Section()
cell.insert('hha2')

stimulus = h.IClamp(0.5)
stimulus.delay = 100
stimulus.dur = 800
stimulus.amp = 7

trec = h.Vector()
vrec = h.Vector()

trec.record(h._ref_t)
vrec.record(cell(0.5)._ref_v)

stimuli_rec = h.Vector()
stimuli_rec.record(stimulus._ref_i)

h.load_file("stdrun.hoc")
h.init()
h.v_init = -65
h.tstop = 1000

before_run = timeit.default_timer()
print('initial time:%f'%(before_run-start))

h.run()

after_run = timeit.default_timer()
print('run time:%f'%(after_run-before_run))

subplot(2,1,1)
plot(trec.as_numpy(),vrec.as_numpy())
xlabel('Time (ms)')
ylabel('Vm (mV)')
title('CA1 Pyramidal Cell')

subplot(2,1,2)
plot(trec.as_numpy(),stimuli_rec.as_numpy())
xlabel('Time (ms)')
ylabel('i (nA)')
title('inject current')

show()




