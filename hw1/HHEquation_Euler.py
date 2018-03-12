# Copyright 2018 dufengtong
#
# a neuron class using HH equation as neuron model
# solving HH equation using Euler method
import matplotlib.pyplot as plt
import numpy as np

class HHNeuron:
    def __init__(self):
        # set parameters
        self.v = -65.0
        self.Cm = 1.0
        self.gl = 0.0
        self.gNa = 120.0
        self.gK = 36.0
        self.ENa = 50.0
        self.EK = -77.0
        self.El = -54.387

        self.dt = 0.05
        self.I = 0.0

        self.t = 0

        self.m = 0
        self.n = 0
        self.h = 0
        self.alpha_m = 0
        self.alpha_n = 0
        self.alpha_h = 0
        self.beta_m = 0
        self.beta_n = 0
        self.beta_h = 0
        self.INa = 0
        self.IK = 0

    def print_params(self):
        print('INa:%f, IK:%f' % (self.INa, self.IK))
        print('m:%f, am:%f, bm:%f' % (self.m, self.alpha_m, self.beta_m))

    def update_all_params(self):
        self.t += self.dt

        self.update_alpha_m()
        self.update_alpha_n()
        self.update_alpha_h()
        self.update_beta_m()
        self.update_beta_n()
        self.update_beta_h()

        self.update_m()
        self.update_n()
        self.update_h()

        self.update_INa()
        self.update_IK()

        # inject current
        self.update_I()

        self.update_v()

    def update_v(self):
        delta_v = (self.gl*(self.El-self.v) + self.INa + self.IK + self.I) / self.Cm
        self.v += delta_v * self.dt

    def update_I(self):
        self.I = 10*(self.t>10) - 10*(self.t>20) + 35*(self.t>30) - 35*(self.t>40)

    def update_INa(self):
        self.INa = self.gNa * (self.m**3) * self.h * (self.ENa-self.v)
    def update_IK(self):
        self.IK = self.gK * (self.n**4) * (self.EK - self.v)

    def update_m(self):
        delta_m = self.alpha_m*(1-self.m) - self.beta_m*self.m
        self.m += delta_m * self.dt
    def update_n(self):
        delta_n = self.alpha_n * (1 - self.n) - self.beta_n * self.n
        self.n += delta_n * self.dt
    def update_h(self):
        delta_h = self.alpha_h * (1 - self.h) - self.beta_h * self.h
        self.h += delta_h * self.dt

    def update_alpha_m(self):
        self.alpha_m = 0.1 * (self.v+40) / \
                       (1 - np.exp(-(self.v+40)/10))
    def update_beta_m(self):
        self.beta_m = 4 * np.exp(-(self.v+65)/18)
    def update_alpha_n(self):
        self.alpha_n = 0.01 * (self.v+55) / \
                       (1 - np.exp(-(self.v+55)/10))
    def update_beta_n(self):
        self.beta_n = 0.125 * np.exp(-(self.v+65)/80)
    def update_alpha_h(self):
        self.alpha_h = 0.07 * np.exp(-(self.v+65)/20)
    def update_beta_h(self):
        self.beta_h = 1 / (np.exp(-(self.v+35)/10) + 1)


if __name__ == '__main__':
    neuron = HHNeuron()

    time = 20 # s
    steps = int(time / neuron.dt) + 1
    x = neuron.dt * np.arange(steps)
    # numerical solution
    y = neuron.v * np.ones(x.shape)

    fig = plt.figure()
    ax = fig.add_subplot(111)


    for i in range(0, steps):
        y[i] = neuron.v
        print('time: %d v: %f ' % (neuron.t, neuron.v))
        neuron.update_all_params()

    ax.set_xlim([0, time])
    ax.set_ylim([-80, 80])
    ax.plot(x, y, color='k')
    plt.show()