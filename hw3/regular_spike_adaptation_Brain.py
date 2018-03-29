# Copyright 2018 Du Fengtong
#
# simulate regular spiking with adaptation

from brian2 import *
import numpy as np

# Parameters
Cm = 1*ufarad
g_na = 100*msiemens
g_k = 80*msiemens
g_m = 10*msiemens
gl = 0.1*msiemens

ENa = 50*mV
EK = -100*mV
El = -67*mV


A = 2.0*uamp

tao_z = 100*ms
eqs = Equations('''
dv/dt = (gl*(El-v)-
         g_na*(m*m*m)*h*(v-ENa)-
         g_k*(n*n*n*n)*(v-EK)-
         g_m*z*(v-EK)+
         I_inj)/Cm : volt
I_inj = A*(t>200*ms)*(t<800*ms) : amp
dm/dt = alpha_m*(1-m)-beta_m*m : 1
dn/dt = alpha_n*(1-n)-beta_n*n : 1
dh/dt = alpha_h*(1-h)-beta_h*h : 1
dz/dt = (z_inf-z)/tao_z : 1
alpha_m = 0.32*(mV**-1)*(54*mV+v)/
         (1.-exp(-(54*mV+v)/(4*mV)))/ms : Hz
beta_m = 0.28*(mV**-1)*(v+27*mV)/
        (exp((v+27*mV)/(5*mV))-1)/ms : Hz
alpha_h = 0.128*exp(-(50*mV+v)/(18*mV))/ms : Hz
beta_h = 4./(1+exp(-(27*mV+v)/(5*mV)))/ms : Hz
alpha_n = 0.032*(mV**-1)*(52*mV+v)/
         (1.-exp(-(52*mV+v)/(5*mV)))/ms : Hz
beta_n = .5*exp(-(57*mV+v)/(40*mV))/ms : Hz
z_inf = 1./(1+exp(-(v+20*mV)/(5*mV))) : 1
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
title('regular spiking with adaptation')
plot(trace.t/ms, trace[0].v/mV)
xlabel('t (ms)')
ylabel('v (mV)')
xlim([0, trace.t[-1]/ms])

# plot ion channel open rate
v = np.arange(-100,100,1)
figure(2)
ax1 = subplot(221)
title('RS with adaptation sodium channel m')
# alpha_m = 0.032*(v+52) / (1 - np.exp(-(v+52)/5))
alpha_m = 0.32*(v+54) / (1 - np.exp(-(v+54)/4))
beta_m = -0.28*(v+27) / (1 - np.exp((v+27)/5))
m_inf = alpha_m / (alpha_m + beta_m)
m_tau = 1 / (alpha_m + beta_m)
ax1.plot(v, m_inf, 'b')
ax1.set_xlabel('v (mV)')
ax1.set_ylabel('inf')
ax2 = ax1.twinx()
ax2.plot(v, m_tau, 'g--')
ax2.set_ylabel('tau(ms)')

ax1 = subplot(222)
title('RS with adaptation sodium channel h')
alpha_h = 0.128*np.exp(-(v+50)/18)
beta_h = 4 / (1 + np.exp(-(v+27)/5))
h_inf = alpha_h / (alpha_h + beta_h)
h_tau = 1 / (alpha_h + beta_h)
ax1.plot(v, h_inf, 'b')
ax1.set_xlabel('v (mV)')
ax1.set_ylabel('inf')
ax2 = ax1.twinx()
ax2.plot(v, h_tau, 'g--')
ax2.set_ylabel('tau(ms)')

ax1 = subplot(223)
title('RS with adaptation potassium channel n')
alpha_n = 0.032*(v+52) / (1 - np.exp(-(v+52)/5))
beta_n = 0.5*np.exp(-(v+57)/40)
n_inf = alpha_n / (alpha_n + beta_n)
n_tau = 1 / (alpha_n + beta_n)
ax1.plot(v, n_inf, 'b')
ax1.set_xlabel('v (mV)')
ax1.set_ylabel('inf')
ax2 = ax1.twinx()
ax2.plot(v, n_tau, 'g--')
ax2.set_ylabel('tau(ms)')

ax1 = subplot(224)
title('RS with adaptation M-type potassium channel z')
z_inf = 1 / (1 + np.exp(-(v+20)/5))
z_tau = 100*np.ones(v.shape)
ax1.plot(v, z_inf, 'b')
ax1.set_xlabel('v (mV)')
ax1.set_ylabel('inf')
ax2 = ax1.twinx()
ax2.plot(v, z_tau, 'g--')
ax2.set_ylabel('tau(ms)')
show()