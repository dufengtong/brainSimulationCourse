# Copyright 2018 Du Fengtong
# test Calcium activated K channel.
# "A computatioNal study on how theta modulated inhibition can account for the long temporal
# windows in the entorhiNal-hippocampal loop". Vassilis C, PaNayiota P. 2015.

from brian2 import *
from brian2.units.constants import faraday_constant as FARADAY, gas_constant as R
import timeit

start = timeit.default_timer()
A = 7000*uamp


'''Soma'''
celsius = 22
Rm = 20000*ohm
gh_soma = 0.00005*siemens
cell = Soma(diameter=10*um)
Cm = 1*ufarad # all compartments use the same Cm
gl = 1/Rm
El = -70*mV
ECa = 140*mV
EK = -80*mV
ENa = 50*mV
Ra = 50

gK = 0.007/5*siemens
gNa = 0.007*siemens
Na_att = 1.0  # Na-conductance attenuation, "s", (ar=1 -> zero attenuation)
tau_m = 0.05*ms
tau_h = 1*ms
tau_n = 3.5*ms
# cell.taumin = 3*ms  # min activation time for "s" attenuation system
gh = gh_soma
Eh = -73*mV  # CA1PC.hoc line156 vhalfl_hd
gA = gh_soma
gM = 0.06*siemens
gCaL = 0.0007*siemens
gCaT = 0.00005*siemens
gsCaR = 0.0003*siemens
gsAHP = 0.0005*siemens
gmAHP = 0.09075*siemens
cac = 0.025*mmolar # middle point of activation fct
cai = 0.0001*mmolar # interNal Ca++ concentration
cao = 2*mmolar # exterNal Ca++ concentration
zetar = 12
gmr = 0.2
eqs = Equations('''
v = -100*mV+0.2*mV*t/ms : volt

# cagk2.mod
# ImAHP = gmAHP*o*(v - EK) : amp
do/dt = (o_inf - o) / tau_o : 1
tau_o = 1 / (a_o + b_o) : second
o_inf = a_o * tau_o : 1
a_o = cai*0.28/(cai + 0.00048*mmolar*exp(-2*0.84*FARADAY/kelvin*v/R/(273.15 + celsius)))/ms : Hz
b_o = 0.48/(1 + cai/(0.00000013*mmolar*exp(-2*1.0*FARADAY/kelvin*v/R/(273.15 + celsius))))/ms : Hz

I_inj = A*(t>100*ms)*(t<400*ms) : amp
''')


P = NeuronGroup(1, model=eqs,
                method='exponential_euler')

# Initialization
# P.v = -65*mV

# Record data
trace = StateMonitor(P, ['v','I_inj','o_inf','tau_o','o'], record=0)
run(1 * second, report='text')

# plot spike train
# figure(1)
# subplot(211)
# title('inject current')
# plot(trace.t/ms, trace[0].I_inj/uA)
# xlabel('t (ms)')
# ylabel('I (uA)')
# xlim([0, trace.t[-1]/ms])
# ylim([0, 2*A/uA])
#
# subplot(212)
# title('regular spiking without adaptation')
# plot(trace.t/ms, trace[0].v/mV)
# xlabel('t (ms)')
# ylabel('v (mV)')
# xlim([0, trace.t[-1]/ms])

figure(1)
subplot(211)
title('o')
plot(trace[0].v/mV, trace[0].o)
xlabel('v')
xlim([-100, 100])


subplot(212)
title('o_inf')
plot(trace[0].v/mV, trace[0].o_inf)
xlabel('v')
xlim([-100, 100])
plt.show()