# Copyright 2018 Du Fengtong
#
# simulate regular spiking without adaptation

from brian2 import *
import numpy as np

# Parameters
Cm = 1*ufarad
g_na = 50*msiemens
g_k = 12*msiemens
g_a = 36*msiemens
g_ca = 2.2*msiemens
gl = 0.05*msiemens

ENa = 35*mV
EK = -75*mV
El = -67*mV
ECa = 120*mV


A = 2.0*uamp

eqs = Equations('''
dv/dt = (gl*(El-v)-g_na*(m*m*m)*h*(v-ENa)-g_k*(n*n*n*n)*(v-EK)- g_a*a*a*a*b*(v-EK)-g_ca*c*c*c*(v-ECa)+I_inj)/Cm : volt
I_inj = A*(t>200*ms)*(t<800*ms) : amp
dm/dt = alpha_m*(1-m)-beta_m*m : 1
dn/dt = alpha_n*(1-n)-beta_n*n : 1
dh/dt = alpha_h*(1-h)-beta_h*h : 1
da/dt = alpha_a*(1-a)-beta_a*a : 1
db/dt = alpha_b*(1-b)-beta_b*b : 1
dc/dt = alpha_c*(1-c)-beta_c*c : 1
alpha_m = 0.6*(mV**-1)*(30*mV+v)/(1.-exp(-(30*mV+v)/(10*mV)))/ms : Hz
beta_m = 20*exp(-(55*mV+v)/(18*mV))/ms : Hz
alpha_h = 0.4*exp(-(50*mV+v)/(20*mV))/ms : Hz
beta_h = 6./(1+exp(-(20*mV+v)/(10*mV)))/ms : Hz
alpha_n = 0.02*(mV**-1)*(40*mV+v)/(1.-exp(-(40*mV+v)/(10*mV)))/ms : Hz
beta_n = 0.4*exp(-(50*mV+v)/(80*mV))/ms : Hz
alpha_a = 0.006*(mV**-1)*(90*mV+v)/(1.-exp(-(90*mV+v)/(10*mV)))/ms : Hz
beta_a = 0.1*exp(-(30*mV+v)/(10*mV))/ms : Hz
alpha_b = 0.04*exp(-(70*mV+v)/(20*mV))/ms : Hz
beta_b = 0.6/(1+exp(-(40*mV+v)/(10*mV)))/ms : Hz
alpha_c = 0.3*(mV**-1)*(13*mV+v)/(1.-exp(-(13*mV+v)/(10*mV)))/ms : Hz
beta_c = 10*exp(-(38*mV+v)/(18*mV))/ms : Hz
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
title('regular spiking without adaptation')
plot(trace.t/ms, trace[0].v/mV)
xlabel('t (ms)')
ylabel('v (mV)')
xlim([0, trace.t[-1]/ms])


# plot ion channel open rate
v = np.arange(-100,100,1)
figure(2)
ax1 = subplot(321)
title('RS without adaptation sodium channel m')
# alpha_m = 0.032*(v+52) / (1 - np.exp(-(v+52)/5))
alpha_m = 0.6*(v+30) / (1 - np.exp(-(v+30)/10))
beta_m = 20*np.exp(-(v+55)/18)
m_inf = alpha_m / (alpha_m + beta_m)
m_tau = 1 / (alpha_m + beta_m)
ax1.plot(v, m_inf, 'b')
ax1.set_xlabel('v (mV)')
ax1.set_ylabel('inf')
ax2 = ax1.twinx()
ax2.plot(v, m_tau, 'g--')
ax2.set_ylabel('tau(ms)')

ax1 = subplot(322)
title('RS without adaptation sodium channel h')
alpha_h = 0.4*np.exp(-(v+50)/20)
beta_h = 6/ (1 + np.exp(-(v+20)/10))
h_inf = alpha_h / (alpha_h + beta_h)
h_tau = 1 / (alpha_h + beta_h)
ax1.plot(v, h_inf, 'b')
ax1.set_xlabel('v (mV)')
ax1.set_ylabel('inf')
ax2 = ax1.twinx()
ax2.plot(v, h_tau, 'g--')
ax2.set_ylabel('tau(ms)')

ax1 = subplot(323)
title('RS without adaptation potassium channel n')
alpha_n = 0.02*(v+40) / (1 - np.exp(-(v+40)/10))
beta_n = 0.4*np.exp(-(v+50)/80)
n_inf = alpha_n / (alpha_n + beta_n)
n_tau = 1 / (alpha_n + beta_n)
ax1.plot(v, n_inf, 'b')
ax1.set_xlabel('v (mV)')
ax1.set_ylabel('inf')
ax2 = ax1.twinx()
ax2.plot(v, n_tau, 'g--')
ax2.set_ylabel('tau(ms)')

ax1 = subplot(324)
title('RS without adaptation A-type potassium channel a')
alpha_a = 0.006*(v+90) / (1 - np.exp(-(v+90)/10))
beta_a = 0.1*np.exp(-(v+30)/10)
a_inf = alpha_a / (alpha_a + beta_a)
a_tau = 1 / (alpha_a + beta_a)
ax1.plot(v, a_inf, 'b')
ax1.set_xlabel('v (mV)')
ax1.set_ylabel('inf')
ax2 = ax1.twinx()
ax2.plot(v, a_tau, 'g--')
ax2.set_ylabel('tau(ms)')

ax1 = subplot(325)
title('RS without adaptation A-type potassium channel b')
alpha_b = 0.04*np.exp(-(v+70)/20)
beta_b = 0.6/(np.exp(-(v+40)/10) +1)
b_inf = alpha_b / (alpha_b + beta_b)
b_tau = 1 / (alpha_b + beta_b)
ax1.plot(v, b_inf, 'b')
ax1.set_xlabel('v (mV)')
ax1.set_ylabel('inf')
ax2 = ax1.twinx()
ax2.plot(v, b_tau, 'g--')
ax2.set_ylabel('tau(ms)')

ax1 = subplot(326)
title('RS without adaptation calcium channel c')
alpha_c = 0.3*(v+13) / (1 - np.exp(-(v+13)/10))
beta_c = 10*np.exp(-(v+38)/18)
c_inf = alpha_c / (alpha_c + beta_c)
c_tau = 1 / (alpha_c + beta_c)
ax1.plot(v, c_inf, 'b')
ax1.set_xlabel('v (mV)')
ax1.set_ylabel('inf')
ax2 = ax1.twinx()
ax2.plot(v, c_tau, 'g--')
ax2.set_ylabel('tau(ms)')
show()