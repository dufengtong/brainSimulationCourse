# Copyright 2018 Du Fengtong
#
# simulate regular spiking with adaptation

from brian2 import *
import numpy as np

# Parameters
Cm = 8.04*ufarad
g_na = 900*msiemens
g_k1 = 1.8*msiemens
g_k3 = 1800*msiemens
gl = 4.1*msiemens

ENa = 60*mV
EK = -90*mV
El = -70*mV


A = 100.0*uamp

eqs = Equations('''
dv/dt = (gl*(El-v)-
         g_na*(m*m*m)*h*(v-ENa)-
         g_k1*(n*n*n*n)*(v-EK)-
         g_k3*p*p*(v-EK)+
         I_inj)/Cm : volt
I_inj = A*(t>200*ms)*(t<800*ms) : amp
dm/dt = alpha_m*(1-m)-beta_m*m : 1
dn/dt = alpha_n*(1-n)-beta_n*n : 1
dh/dt = alpha_h*(1-h)-beta_h*h : 1
dp/dt = alpha_p*(1-p)-beta_p*p : 1
alpha_m = 40*(mV**-1)*(v-75*mV)/
         (1.-exp(-(v-75*mV)/(13.5*mV)))/ms : Hz
beta_m = 1.2262*exp(-v/(42.248*mV))/ms : Hz
alpha_h = 0.0035*exp(-v/(24.186*mV))/ms : Hz
beta_h = 0.017*(mV**-1)*(v+51.25*mV)/
         (1.-exp(-(v+51.25*mV)/(5.2*mV)))/ms : Hz
alpha_n = 0.014*(mV**-1)*(44*mV+v)/
         (1.-exp(-(44*mV+v)/(2.3*mV)))/ms : Hz
beta_n = 0.0043*exp(-(44*mV+v)/(34*mV))/ms : Hz
alpha_p = 1*(mV**-1)*(v-95*mV)/
         (1.-exp(-(v-95*mV)/(11.8*mV)))/ms : Hz
beta_p = 0.025*exp(-v/(22.222*mV))/ms : Hz
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
title('fast spiking with adaptation')
plot(trace.t/ms, trace[0].v/mV)
xlabel('t (ms)')
ylabel('v (mV)')
xlim([0, trace.t[-1]/ms])

# plot ion channel open rate
v = np.arange(-100,100,1)
figure(2)
ax1 = subplot(221)
title('FS with adaptation sodium channel m')
# alpha_m = 0.032*(v+52) / (1 - np.exp(-(v+52)/5))
alpha_m = 40*(v-75) / (1 - np.exp(-(v-75)/13.5))
beta_m = 1.2262*np.exp(-v/42.248)
m_inf = alpha_m / (alpha_m + beta_m)
m_tau = 1 / (alpha_m + beta_m)
ax1.plot(v, m_inf, 'b')
ax1.set_xlabel('v (mV)')
ax1.set_ylabel('inf')
ax2 = ax1.twinx()
ax2.plot(v, m_tau, 'g--')
ax2.set_ylabel('tau(ms)')

ax1 = subplot(222)
title('FS with adaptation sodium channel h')
alpha_h = 0.0035*np.exp(-v/24.186)
beta_h = 0.017*(v+51.25) / (1 - np.exp(-(v+51.25)/5.2))
h_inf = alpha_h / (alpha_h + beta_h)
h_tau = 1 / (alpha_h + beta_h)
ax1.plot(v, h_inf, 'b')
ax1.set_xlabel('v (mV)')
ax1.set_ylabel('inf')
ax2 = ax1.twinx()
ax2.plot(v, h_tau, 'g--')
ax2.set_ylabel('tau(ms)')

ax1 = subplot(223)
title('FS with adaptation k1 potassium channel n')
alpha_n = 0.014*(v+44) / (1 - np.exp(-(v+44)/2.3))
beta_n = 0.0043*np.exp(-(v+44)/34)
n_inf = alpha_n / (alpha_n + beta_n)
n_tau = 1 / (alpha_n + beta_n)
ax1.plot(v, n_inf, 'b')
ax1.set_xlabel('v (mV)')
ax1.set_ylabel('inf')
ax2 = ax1.twinx()
ax2.plot(v, n_tau, 'g--')
ax2.set_ylabel('tau(ms)')

ax1 = subplot(224)
title('FS with adaptation k3 potassium channel p')
alpha_p = 1*(v-95) / (1 - np.exp(-(v-95)/11.8))
beta_p = 0.025*np.exp(-v/22.222)
p_inf = alpha_p / (alpha_p + beta_p)
p_tau = 1 / (alpha_p + beta_p)
ax1.plot(v, p_inf, 'b')
ax1.set_xlabel('v (mV)')
ax1.set_ylabel('inf')
ax2 = ax1.twinx()
ax2.plot(v, p_tau, 'g--')
ax2.set_ylabel('tau(ms)')
show()
show()