# Copyright 2018 Du Fengtong
#
# simulate regular spiking without adaptation

from brian2 import *
import numpy as np

# Parameters
Cm = 1*ufarad
g_na = 35*msiemens
g_k = 6*msiemens
g_a = 1.4*msiemens
g_nap = 0.2*msiemens
g_m = 3*msiemens
gl = 0.05*msiemens

ENa = 55*mV
EK = -90*mV
El = -70*mV

b_tau = 15*ms
z_tau = 75*ms

A = 1.7*uamp

eqs = Equations('''
dv/dt = (gl*(El-v)-
        g_na*(m_inf*m_inf*m_inf)*h*(v-ENa)-
        g_k*(n*n*n*n)*(v-EK)-
        g_a*(a_inf*a_inf*a_inf)*b*(v-EK)-
        g_nap*p_inf*(v-ENa)-
        g_m*z*(v-EK)+
        I_inj)/Cm : volt
I_inj = A*(t>200*ms)*(t<800*ms) : amp
dn/dt = (n_inf-n) / n_tau : 1
dh/dt = (h_inf-h) / h_tau : 1
db/dt = (b_inf-b) / b_tau : 1
dz/dt = (z_inf-z) / z_tau : 1
m_inf = 1./(1+exp(-(30*mV+v)/(9.5*mV))) : 1
h_inf = 1./(1+exp((45*mV+v)/(7*mV))) : 1
h_tau = (0.1 + 0.75/(1+exp((40.5*mV+v)/(6*mV))))*ms : second
n_inf = 1./(1+exp(-(35*mV+v)/(10*mV))) : 1
n_tau = (0.1 + 0.5/(1+exp((27*mV+v)/(15*mV))))*ms : second
p_inf = 1./(1+exp(-(45*mV+v)/(3*mV))) : 1
a_inf = 1./(1+exp(-(50*mV+v)/(20*mV))) : 1
b_inf = 1./(1+exp((80*mV+v)/(6*mV))) : 1
z_inf = 1./(1+exp(-(39*mV+v)/(5*mV))) : 1
''')

P = NeuronGroup(1, model=eqs, threshold='v>-20*mV',
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
title('bursting excitatory neuron')
plot(trace.t/ms, trace[0].v/mV)
xlabel('t (ms)')
ylabel('v (mV)')
xlim([0, trace.t[-1]/ms])


# plot ion channel open rate
# v = np.arange(-100,100,1)
# figure(2)
# ax1 = subplot(321)
# title('RS without adaptation sodium channel m')
# # alpha_m = 0.032*(v+52) / (1 - np.exp(-(v+52)/5))
# alpha_m = 0.6*(v+30) / (1 - np.exp(-(v+30)/10))
# beta_m = 20*np.exp(-(v+55)/18)
# m_inf = alpha_m / (alpha_m + beta_m)
# m_tau = 1 / (alpha_m + beta_m)
# ax1.plot(v, m_inf, 'b')
# ax1.set_xlabel('v (mV)')
# ax1.set_ylabel('inf')
# ax2 = ax1.twinx()
# ax2.plot(v, m_tau, 'g--')
# ax2.set_ylabel('tau(ms)')
#
# ax1 = subplot(322)
# title('RS without adaptation sodium channel h')
# alpha_h = 0.4*np.exp(-(v+50)/20)
# beta_h = 6/ (1 + np.exp(-(v+20)/10))
# h_inf = alpha_h / (alpha_h + beta_h)
# h_tau = 1 / (alpha_h + beta_h)
# ax1.plot(v, h_inf, 'b')
# ax1.set_xlabel('v (mV)')
# ax1.set_ylabel('inf')
# ax2 = ax1.twinx()
# ax2.plot(v, h_tau, 'g--')
# ax2.set_ylabel('tau(ms)')
#
# ax1 = subplot(323)
# title('RS without adaptation potassium channel n')
# alpha_n = 0.02*(v+40) / (1 - np.exp(-(v+40)/10))
# beta_n = 0.4*np.exp(-(v+50)/80)
# n_inf = alpha_n / (alpha_n + beta_n)
# n_tau = 1 / (alpha_n + beta_n)
# ax1.plot(v, n_inf, 'b')
# ax1.set_xlabel('v (mV)')
# ax1.set_ylabel('inf')
# ax2 = ax1.twinx()
# ax2.plot(v, n_tau, 'g--')
# ax2.set_ylabel('tau(ms)')
#
# ax1 = subplot(324)
# title('RS without adaptation A-type potassium channel a')
# alpha_a = 0.006*(v+90) / (1 - np.exp(-(v+90)/10))
# beta_a = 0.1*np.exp(-(v+30)/10)
# a_inf = alpha_a / (alpha_a + beta_a)
# a_tau = 1 / (alpha_a + beta_a)
# ax1.plot(v, a_inf, 'b')
# ax1.set_xlabel('v (mV)')
# ax1.set_ylabel('inf')
# ax2 = ax1.twinx()
# ax2.plot(v, a_tau, 'g--')
# ax2.set_ylabel('tau(ms)')
#
# ax1 = subplot(325)
# title('RS without adaptation A-type potassium channel b')
# alpha_b = 0.04*np.exp(-(v+70)/20)
# beta_b = 0.6/(np.exp(-(v+40)/10) +1)
# b_inf = alpha_b / (alpha_b + beta_b)
# b_tau = 1 / (alpha_b + beta_b)
# ax1.plot(v, b_inf, 'b')
# ax1.set_xlabel('v (mV)')
# ax1.set_ylabel('inf')
# ax2 = ax1.twinx()
# ax2.plot(v, b_tau, 'g--')
# ax2.set_ylabel('tau(ms)')
#
# ax1 = subplot(326)
# title('RS without adaptation calcium channel c')
# alpha_c = 0.3*(v+13) / (1 - np.exp(-(v+13)/10))
# beta_c = 10*np.exp(-(v+38)/18)
# c_inf = alpha_c / (alpha_c + beta_c)
# c_tau = 1 / (alpha_c + beta_c)
# ax1.plot(v, c_inf, 'b')
# ax1.set_xlabel('v (mV)')
# ax1.set_ylabel('inf')
# ax2 = ax1.twinx()
# ax2.plot(v, c_tau, 'g--')
# ax2.set_ylabel('tau(ms)')
show()