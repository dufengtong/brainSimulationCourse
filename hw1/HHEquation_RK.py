import matplotlib.pyplot as plt
import numpy as np
from HHEquation_Euler import HHNeuron

def runge_kutta(f, y, dt):
    y1 = f(y) * dt
    y2 = f(y + y1 / 2) * dt
    y3 = f(y + y2 / 2) * dt
    y4 = f(y + y3) * dt
    return y + (y1 + 2 * y2 + 2 * y3 + y4) / 6

class HHNeuron_RK:
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
        self.v = runge_kutta(self.d_v, self.v, self.dt)
    def update_I(self):
        self.I = 10 * (self.t > 10) - 10 * (self.t > 20) + 35 * (self.t > 30) - 35 * (self.t > 40)
    def update_INa(self):
        self.INa = self.gNa * (self.m ** 3) * self.h * (self.ENa - self.v)
    def update_IK(self):
        self.IK = self.gK * (self.n ** 4) * (self.EK - self.v)
    def update_m(self):
        self.m = runge_kutta(self.d_m, self.m, self.dt)
    def update_n(self):
        self.n = runge_kutta(self.d_n, self.n, self.dt)
    def update_h(self):
        self.h = runge_kutta(self.d_h, self.h, self.dt)
    def update_alpha_m(self):
        self.alpha_m = 0.1 * (self.v + 40) / \
        (1 - np.exp(-(self.v + 40) / 10))
    def update_beta_m(self):
        self.beta_m = 4 * np.exp(-(self.v + 65) / 18)
    def update_alpha_n(self):
        self.alpha_n = 0.01 * (self.v + 55) / \
        (1 - np.exp(-(self.v + 55) / 10))
    def update_beta_n(self):
        self.beta_n = 0.125 * np.exp(-(self.v + 65) / 80)
    def update_alpha_h(self):
        self.alpha_h = 0.07 * np.exp(-(self.v + 65) / 20)
    def update_beta_h(self):
        self.beta_h = 1 / (np.exp(-(self.v + 35) / 10) + 1)
    def d_v(self, v):
        return (self.gl * (self.El - v) + self.INa + self.IK + self.I) / self.Cm
    def d_m(self, m):
        return self.alpha_m * (1 - m) - self.beta_m * m
    def d_n(self, n):
        return self.alpha_n * (1 - n) - self.beta_n * n
    def d_h(self, h):
        return self.alpha_h * (1 - h) - self.beta_h * h

if __name__ == '__main__':
    rk_neuron = HHNeuron_RK()
    euler_neuron = HHNeuron()
    time = 20 # s
    steps = int(time / rk_neuron.dt) + 1
    x = rk_neuron.dt * np.arange(steps)
    # numerical solution
    rk_y = rk_neuron.v * np.ones(x.shape)
    euler_y = euler_neuron.v * np.ones(x.shape)
    fig = plt.figure()
    ax = fig.add_subplot(111)

    for i in range(0, steps):
        rk_y[i] = rk_neuron.v
        euler_y[i] = euler_neuron.v
        print('time: %d rk_v: %f euler_v: %f diff: %f' % (rk_neuron.t, rk_neuron.v, euler_neuron.v,
        rk_neuron.v-euler_neuron.v))
        rk_neuron.update_all_params()
        euler_neuron.update_all_params()

    ax.set_xlim([0, time])
    ax.set_ylim([-80, 80])
    rk_handle, = ax.plot(x, rk_y, color='r', label='runge kutta')
    euler_handle, = ax.plot(x, euler_y, color='y', label='euler')
    plt.legend(handles=[euler_handle, rk_handle])
    plt.show()