# Copyright 2018 Du Fengtong
#
# simulate regular spiking with adaptation

from brian2 import *
import numpy as np

# Parameters
Cm = 1*ufarad
g_na = 35*msiemens
g_k = 9*msiemens
gl = 0.1*msiemens

ENa = 55*mV
EK = -90*mV
El = -65*mV


A = 2.0*uamp

eqs = Equations('''
dv/dt = (gl*(El-v)-
         g_na*(m_inf*m_inf*m_inf)*h*(v-ENa)-
         g_k*(n*n*n*n)*(v-EK)+
         I_inj)/Cm : volt
I_inj = A*(t>200*ms)*(t<800*ms) : amp
dn/dt = alpha_n*(1-n)-beta_n*n : 1
dh/dt = alpha_h*(1-h)-beta_h*h : 1
m_inf = alpha_m / (alpha_m+beta_m) : 1
alpha_m = 0.1*(mV**-1)*(35*mV+v)/
         (1.-exp(-(35*mV+v)/(10*mV)))/ms : Hz
beta_m = 4*exp(-(60*mV+v)/(18*mV))/ms : Hz
alpha_h = 0.07*exp(-(58*mV+v)/(20*mV))/ms : Hz
beta_h = 1./(1+exp(-(28*mV+v)/(10*mV)))/ms : Hz
alpha_n = 0.01*(mV**-1)*(34*mV+v)/
         (1.-exp(-(34*mV+v)/(10*mV)))/ms : Hz
beta_n = 0.125*exp(-(44*mV+v)/(80*mV))/ms : Hz
''')

P = NeuronGroup(1, model=eqs, threshold='v>-20*mV', refractory=3*ms,
                method='exponential_euler')

# Initialization
P.v = 'El + (randn() * 5 - 5)*mV'

# Record data
trace = StateMonitor(P, ['v','I_inj'], record=0)
run(1 * second, report='text')

# plot spike train
figure(1)
subplot(211)
title('inject current')
plot(trace.t/ms, trace[0].I_inj/uA)
xlabel('t (ms)')
ylabel('I (uA)')
xlim([0, trace.t[-1]/ms])
ylim([0, 2*A/uA])

subplot(212)
title('fast spiking without adaptation')
plot(trace.t/ms, trace[0].v/mV)
xlabel('t (ms)')
ylabel('v (mV)')
xlim([0, trace.t[-1]/ms])

# plot ion channel open rate
v = np.arange(-100,100,1)
figure(2)
ax1 = subplot(131)
title('FS without adaptation sodium channel m')
alpha_m = 0.1*(v+35) / (1 - np.exp(-(v+35)/10))
beta_m = 4*np.exp(-(v+60)/18)
m_inf = alpha_m / (alpha_m + beta_m)
m_tau = 1 / (alpha_m + beta_m)
ax1.plot(v, m_inf, 'b')
ax1.set_xlabel('v (mV)')
ax1.set_ylabel('inf')
ax2 = ax1.twinx()
ax2.plot(v, m_tau, 'g--')
ax2.set_ylabel('tau(ms)')

ax1 = subplot(132)
title('FS without adaptation sodium channel h')
alpha_h = 0.07*np.exp(-(v+58)/20)
beta_h = 1 / (1 + np.exp(-(v+28)/10))
h_inf = alpha_h / (alpha_h + beta_h)
h_tau = 1 / (alpha_h + beta_h)
ax1.plot(v, h_inf, 'b')
ax1.set_xlabel('v (mV)')
ax1.set_ylabel('inf')
ax2 = ax1.twinx()
ax2.plot(v, h_tau, 'g--')
ax2.set_ylabel('tau(ms)')

ax1 = subplot(133)
title('FS without adaptation potassium channel n')
alpha_n = 0.01*(v+34) / (1 - np.exp(-(v+34)/10))
beta_n = 0.125*np.exp(-(v+44)/80)
n_inf = alpha_n / (alpha_n + beta_n)
n_tau = 1 / (alpha_n + beta_n)
ax1.plot(v, n_inf, 'b')
ax1.set_xlabel('v (mV)')
ax1.set_ylabel('inf')
ax2 = ax1.twinx()
ax2.plot(v, n_tau, 'g--')
ax2.set_ylabel('tau(ms)')

show()