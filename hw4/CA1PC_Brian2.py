# Copyright 2018 Du Fengtong
# simulate CA1 Pyramidal Cell based on
# "A computatioNal study on how theta modulated inhibition can account for the long temporal
# windows in the entorhiNal-hippocampal loop". Vassilis C, PaNayiota P. 2015.

from brian2 import *
import numpy as np
from brian2.units.constants import faraday_constant as FARADAY, gas_constant as R
import timeit
from mpl_toolkits.mplot3d import Axes3D
#
# cell.axon = Cylinder(length=150*um, diameter=1*um, n=1)
#
# cell.dendrite = Cylinder(length=50*um, diameter=2*um, n=5)
start = timeit.default_timer()
A = 60*uA  # 70*mamp  # mA


'''Soma'''
celsius = 22
Rm = 20000*ohm
gh_soma = 0.00005*siemens
cell = Soma(diameter=10*um)
# soma_surface = (4/3)*pi*(10/2)**3 # 10um=0.01mm=0.0001cm
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
gmAHP = 0.01*siemens
cac = 0.025*mmolar # middle point of activation fct
# cai = 0.0001*mmolar # inf Ca++ concentration
cao = 2*mmolar # exterNal Ca++ concentration
zetar = 12
gmr = 0.2
eqs = Equations('''
dv/dt = (-IL - INa - IK - IA- IM - ICaL
        - Ibuff
        - ICaT
        - IsCaR
        - IsAHP
        - ImAHP
        + I_inj)/Cm : volt
# hha2.mod
IL = gl*(v - El) : amp
INa = gNa*(m*m)*h*s*(v - ENa) : amp
dm/dt = (m_inf - m)*(1+exp(-dt/tau_m))/dt : 1
m_inf = 1./(1.+exp(-(44.*mV+v)/(3*mV))) : 1
dh/dt = (1-exp(-dt/tau_h))*(h_inf - h)/dt : 1
h_inf = 1 / (1+exp((v+49*mV)/(3.5*mV))) : 1
ds/dt = (1-exp(-dt/tau_s))*(s_inf - s)/dt : 1
s_inf = (1 + exp((v+73*mV)/(2*mV))) / (1+exp((v+73*mV)/(2*mV))) : 1
tau_s_tmp = exp(0.001*zetar*gmr*(mV**-1)*(v+73*mV)*96480/(8.315*(273.16+celsius))) /
            (1+exp(0.001*zetar*(mV**-1)*(v+73*mV)*96480/(8.315*(273.16+celsius))))/(0.0003/ms)  : second
tau_s = (tau_s_tmp<3*ms)*(tau_s_tmp - 3*ms) + 3*ms : second

IK = gK*(n*n)*(v - EK) : amp
dn/dt = (1-exp(-dt/tau_n))*(n_inf - n)/dt : 1
n_inf = 1 / (1+exp(-(v+46.3*mV)/(3*mV))) : 1

# kap.mod
IA = gA*(n_A*n_A*n_A*n_A)*l*(v-EK) : amp
dn_A/dt = (1 - exp(-dt/(0.2*ms)))*(n_A_inf - n_A)/dt : 1
n_A_inf = a_n_A /(a_n_A + b_n_A) : 1
tau_l = (v>(20*mV))*2.6*ms*(v+20*mV)/(10*mV) + 5*ms : second
a_n_A = -0.01*(mV**-1)*(v+21.3*mV)/(exp((v+21.3*mV)/(-35*mV))-1)/ms : Hz
b_n_A = 0.01*(mV**-1)*(v+21.3*mV)/(exp((v+21.3*mV)/(35*mV))-1)/ms : Hz
dl/dt = (1 - exp(-dt/tau_l))*(l_inf - l)/dt : 1
l_inf = a_l /(a_l + b_l) : 1
a_l = -0.01*(mV**-1)*(v+58*mV)/(exp((v+58*mV)/(8.2*mV))-1)/ms : Hz
b_l = 0.01*(mV**-1)*(v+58*mV)/(exp((v+58*mV)/(-8.2*mV))-1)/ms : Hz

# km.mod
IM = 0.0001*2.3**((celsius - 23)/10)*gM*n_M*(v - EK) : amp
dn_M/dt = (1 - exp(-dt*2.3**((celsius-23)/10)/tau_n_M))*(n_M_inf-n_M)/dt : 1
n_M_inf = a_n_M * tau_n_M : 1
tau_n_M = 1 / (a_n_M + b_n_M) : second
a_n_M = 0.001*(mV**-1)*(v+30*mV) / (1 - exp(-(v+30*mV)/(9*mV)))/ms : Hz
b_n_M = -0.001*(mV**-1)*(v+30*mV) / (1 - exp((v+30*mV)/(9*mV)))/ms : Hz

# cal.mod
ICaL = gCaL*mCaL*(0.001*mmolar/(0.001*mmolar+cai))*ghk : amp
dmCaL/dt = (1 - exp(-dt/tau_mCaL))*(mCaL_inf - mCaL)/dt : 1
mCaL_inf = a_mCaL /(a_mCaL + b_mCaL) : 1
tau_mCaL = 1 / (5 * (a_mCaL + b_mCaL)) : second
a_mCaL = 0.055*(mV**-1)*(-27.01*mV - v)/(exp((-27.01*mV-v)/(3.8*mV)) - 1)/ms : Hz
b_mCaL = 0.94*exp((-63.01*mV-v)/(17*mV))/ms : Hz
KTF = ((25*mV/293.15)*(celsius + 273.15)) : volt
nu = 2*v / KTF : 1
efun = (abs(nu)<0.0001)*(1-nu/2) + (abs(nu)>0.0001)*(nu/(exp(nu) - 1)) : 1
ghk = (-(KTF/2)*(1. - (cai/cao)*exp(nu))*efun) : volt

# cat.mod
ICaT = gCaT*mCaT*mCaT*hCaT*(0.001*mmolar/(0.001*mmolar+cai))*ghk : amp
dmCaT/dt = (1 - exp(-dt/tau_mCaT))*(mCaT_inf - mCaT)/dt : 1
mCaT_inf = a_mCaT /(a_mCaT + b_mCaT) : 1
tau_mCaT = 1 / (a_mCaT + b_mCaT) : second
a_mCaT = 0.1967*(mV**-1)*(-v+19.88*mV)/(exp((-v+19.88*mV)/(10.0*mV))-1.0)/ms : Hz
b_mCaT = 0.046*exp(-v/(22.73*mV))/ms : Hz
dhCaT/dt = (1 - exp(-dt/tau_hCaT))*(hCaT_inf - hCaT)/dt : 1
hCaT_inf = a_hCaT /(a_hCaT + b_hCaT) : 1
tau_hCaT = 1 / (0.68 * (a_hCaT + b_hCaT)) : second
a_hCaT = 0.00016*exp(-(v+57*mV)/(19*mV))/ms : Hz
b_hCaT = 1/(exp((-v+15*mV)/(10*mV))+1.0)/ms : Hz

# cad.mod
Ibuff = - ICaL
        - ICaT
        - IsCaR
        : amp
dcai/dt = drive_channel*mmolar/18/ms + (0.0001*mmolar -cai)/(200*ms)*7 : mmolar
drive_channel_tmp =  - (10000) * Ibuff / (2 * FARADAY*katal * 0.00001)  : 1
drive_channel = (drive_channel_tmp>0)*drive_channel_tmp : 1


# somacar.mod
IsCaR = gsCaR*msCaR*msCaR*msCaR*hsCaR*(v - ECa) : amp
dmsCaR/dt = (1 - exp(-dt/(100*ms)))*(msCaR_inf - msCaR)/dt : 1
dhsCaR/dt = (1 - exp(-dt/(5*ms)))*(hsCaR_inf - hsCaR)/dt : 1
msCaR_inf = 1 / (1 + exp((v+60*mV)/(-3*mV))) : 1
hsCaR_inf = 1/ (1 + exp((v+62*mV)/(1*mV))) : 1

# kca.mod
IsAHP = gsAHP*msAHP*msAHP*msAHP*(v - EK) : amp
dmsAHP/dt = (msAHP_inf - msAHP) / tau_msAHP : 1
msAHP_inf = car / ( 1 + car ) : 1
tau_msAHP_tmp = 1 / (0.03/ms) / (1 + car) / (3**((celsius-22.0)/10)) : second
tau_msAHP = (tau_msAHP_tmp>0.5*ms)*(tau_msAHP_tmp - 0.5*ms) + 0.5*ms : second
car = (cai/cac)**2 : 1

# cagk2.mod
ImAHP = gmAHP*o*(v - EK) : amp
do/dt = (o_inf - o) / tau_o : 1
tau_o = 1 / (a_o + b_o) : second
o_inf = a_o * tau_o : 1
a_o = cai*0.28/(cai + 0.00048*mmolar*exp(-2*0.84*FARADAY/kelvin*v/R/(273.15 + celsius)))/ms : Hz
b_o = 0.48/(1 + cai/(0.00000013*mmolar*exp(-2*1.0*FARADAY/kelvin*v/R/(273.15 + celsius))))/ms : Hz

# I_inj = A*(t>200*ms)*(t<800*ms)/soma_surface : amp
I_inj = A*(t>200*ms)*(t<800*ms) : amp
''')


P = NeuronGroup(1, model=eqs,
                method='exponential_euler')

# Initialization
P.v = -65*mV

# Record data
trace = StateMonitor(P, ['v','I_inj'], record=0)

before_run = timeit.default_timer()
print('initial time:%f'%(before_run-start))

run(1 * second, report='text')

after_run = timeit.default_timer()
print('run time:%f'%(after_run-before_run))

# plot spike train
figure(1)
subplot(211)
title('inject current')
plot(trace.t/ms, trace[0].I_inj/uA)
xlabel('t (ms)')
ylabel('I (uA)')
xlim([0, trace.t[-1]/ms])

subplot(212)
title('CA1 Pyramidal Soma')
plot(trace.t/ms, trace[0].v/mV)
xlabel('t (ms)')
ylabel('v (mV)')


# figure(2)
# subplot(211)
# title('o')
# plot(trace[0].v/mV, trace[0].o)
# xlabel('v')
# xlim([-100, 100])
#
#
# subplot(212)
# title('o_inf')
# plot(trace[0].v/mV, trace[0].o_inf)
# xlabel('v')
# xlim([-100, 100])

plt.show()