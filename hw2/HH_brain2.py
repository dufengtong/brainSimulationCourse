# Copyright 2018 Du Fengtong
#
# This code is written based on Brain2 official example "COBAHH"
# parameters are following paper 'A review of tools and strategies' (2006). Brette

from brian2 import *

# Parameters
Cm = 1*ufarad
gl = 5e-5*siemens

El = -60*mV
EK = -90*mV
ENa = 50*mV
g_na = 100*msiemens
g_kd = 30*msiemens
VT = -63*mV
A = 0.4*uamp
f = 10*Hz

# The inject current equations:
# direct current
# I_inj = A : amp
# sin current
# I_inj = A*(sin(2*pi*f*t)+1) : amp
# ramp current
# I_inj = A/second*t : amp
# noise current
# I_inj = A*(1 + 0.1*randn()) : amp (constant over dt)
tao = 1000*ms
eqs = Equations('''
dv/dt = (gl*(El-v)-
         g_na*(m*m*m)*h*(v-ENa)-
         g_kd*(n*n*n*n)*(v-EK)+
         I_inj)/Cm : volt
I_inj = A*(1 + 0.1*randn()) : amp (constant over dt)
dm/dt = alpha_m*(1-m)-beta_m*m : 1
dn/dt = alpha_n*(1-n)-beta_n*n : 1
dh/dt = alpha_h*(1-h)-beta_h*h : 1
alpha_m = 0.32*(mV**-1)*(13*mV-v+VT)/
         (exp((13*mV-v+VT)/(4*mV))-1.)/ms : Hz
beta_m = 0.28*(mV**-1)*(v-VT-40*mV)/
        (exp((v-VT-40*mV)/(5*mV))-1)/ms : Hz
alpha_h = 0.128*exp((17*mV-v+VT)/(18*mV))/ms : Hz
beta_h = 4./(1+exp((40*mV-v+VT)/(5*mV)))/ms : Hz
alpha_n = 0.032*(mV**-1)*(15*mV-v+VT)/
         (exp((15*mV-v+VT)/(5*mV))-1.)/ms : Hz
beta_n = .5*exp((10*mV-v+VT)/(40*mV))/ms : Hz
''')

P = NeuronGroup(1, model=eqs, threshold='v>-20*mV', refractory=3*ms,
                method='exponential_euler')

# Initialization
P.v = 'El + (randn() * 5 - 5)*mV'

# Record data
trace = StateMonitor(P, ['v','I_inj'], record=0)
run(1 * second, report='text')

# plot results
subplot(211)
title('inject current')
plot(trace.t/ms, trace[0].I_inj/uA)
xlabel('t (ms)')
ylabel('I (uA)')
xlim([0, trace.t[-1]/ms])
ylim([0, 2*A/uA])

subplot(212)
title('HH neuron with exponential euler method')
plot(trace.t/ms, trace[0].v/mV)
xlabel('t (ms)')
ylabel('v (mV)')
xlim([0, trace.t[-1]/ms])
show()